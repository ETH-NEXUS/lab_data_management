import asyncio
import logging

import websockets
from django.conf import settings

from channels.exceptions import DenyConnection
from channels.generic.websocket import AsyncWebsocketConsumer
from django.utils.functional import cached_property

logger = logging.getLogger(__name__)


class WebsocketProxyConsumer(AsyncWebsocketConsumer):
    """
    Proxying websocket connections.
    Found here: https://gist.github.com/brianglass/e3184341afe63ed348144753ee62dce5
    """

    # This is the frequency of pinging we do to the target url.  Pinging seems
    # to confuse the code-server connection and it loses connection every 20
    # seconds, so for now we'll default to no pinging.
    PING_INTERVAL = None

    # This is the maximum size of frames going to/from the target url. We have
    # seen some frames larger than 1MiB being sent between the VS Code client and
    # code-server.
    MAX_SIZE = 2097152  # 2 MiB

    # These headers are passed through from the client to the target url.
    PASSTHROUGH_HEADERS = {}

    async def get_target_url(self):
        """This should be overriden in child classes."""
        ws_url = settings.JUPYTER_URL.replace("https://", "wss://").replace(
            "http://", "ws://"
        )
        uri = self.scope["url_route"]["kwargs"]["uri"]
        query_string = self.scope["query_string"].decode("utf-8")
        url = f"{ws_url}/{uri}?{query_string}"
        return url

    async def connect(self):
        """Establish connections to both the client and the target url."""

        target_url = await self.get_target_url()

        # The requested url is not valid.
        if target_url is None:
            logger.warning("WebSocketProxy Denying websocket connection.")
            raise DenyConnection("The requested endpoint is not valid.")

        # Connect to the target url.
        try:
            logger.debug(f"WebSocketProxy Connecting to {target_url}...")
            # This dows the trick if you have security enabled on the Jupyter notebook
            headers = self.passthrough_headers
            headers.update({"Authorization": f"token {settings.NOTEBOOK_SECURE_TOKEN}"})
            self.websocket = await websockets.connect(
                target_url,
                max_size=self.MAX_SIZE,
                ping_interval=self.PING_INTERVAL,
                extra_headers=self.passthrough_headers,
                subprotocols=self.scope["subprotocols"],
                origin=self.request_headers.get("Origin"),
            )
        except websockets.InvalidURI:
            logger.exception(
                "WebSocketProxy The requested endpoint could not be reached."
            )
            raise DenyConnection(
                "WebSocketProxy The requested endpoint could not be reached."
            )
        except websockets.InvalidHandshake:
            logger.exception(
                "WebSocketProxy Communication with the target url was incoherent."
            )
            raise DenyConnection(
                "WebSocketProxy Communication with the target url was incoherent."
            )

        # Accept the client connection. Use the subprotocol negotiated with the
        # target url.
        await self.accept(self.websocket.subprotocol)

        # Forward packets from the target websocket back to the client.
        self.consumer_task = asyncio.create_task(self.consume_from_target())

    @cached_property
    def request_headers(self):
        return {
            h.decode("utf-8").title(): v.decode("utf-8")
            for h, v in self.scope["headers"]
        }

    @cached_property
    def passthrough_headers(self):
        return {
            h: v
            for h, v in self.request_headers.items()
            if h in self.PASSTHROUGH_HEADERS
        }

    async def disconnect(self, close_code):
        """The websocket consumer is shutting down. Shut down the connection to
        the target url."""

        # Disconnect can be called before self.consumer_task is created.

        if hasattr(self, "consumer_task"):
            self.consumer_task.cancel()

            # Let the task complete
            await self.consumer_task

    async def receive(self, text_data=None, bytes_data=None):
        """Forward packets from the client to the target url."""

        try:
            await self.websocket.send(bytes_data or text_data)
        except websockets.ConnectionClosedError:
            # The target probably closed the connection.
            logger.exception(
                "WebSocketProxy The outgoing connection was closed by the target."
            )
            await self.close()

    async def consume_from_target(self):
        """A websocket consumer to forward data from the target url to the client."""

        try:
            async for data in self.websocket:
                if hasattr(data, "decode"):
                    await self.send(bytes_data=data)
                else:
                    await self.send(text_data=data)
        except asyncio.exceptions.CancelledError:
            # This is triggered by the consumer itself when the client connection is terminating.
            logger.debug(
                "Shutting down the websocket consumer task and closing the outgoing websocket."
            )
            await self.websocket.close()
        except websockets.ConnectionClosedError:
            # The target probably closed the connection.
            logger.exception(
                "WebSocketProxy The outgoing connection was closed by the target."
            )
            await self.close()
