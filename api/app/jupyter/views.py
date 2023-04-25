from django.conf import settings
from revproxy.views import ProxyView
from django.core.handlers.wsgi import WSGIRequest
from django.http import HttpResponseRedirect


# We need to replace quote_plus with quote to replace forwarded utls with %20 instead of +
# Jupyter understands only %20
from urllib.parse import quote as quote_plus

# Chars that don't need to be quoted. We use same than nginx:
#   https://github.com/nginx/nginx/blob/nginx-1.9/src/core/ngx_string.c
#   (Lines 1433-1449)
QUOTE_SAFE = r'<.;>\(}*+|~=-$/_:^@)[{]&\'!,"`'


class JupyterProxyView(ProxyView):
    upstream = settings.JUPYTER_URL

    def get_quoted_path(self, path):
        """Return quoted path to be used in proxied request"""
        # We need to replace quote_plus with quote to replace forwarded utls with %20 instead of +
        # Jupyter understands only %20
        return quote_plus(path.encode("utf8"), QUOTE_SAFE)

    def get_request_headers(self):
        request_headers = super().get_request_headers()
        request_headers.update(
            {"Authorization": f"token {settings.NOTEBOOK_SECURE_TOKEN}"}
        )
        return request_headers

    def dispatch(self, request: WSGIRequest, path):
        if request.user.is_authenticated:
            return super().dispatch(request, path)
        else:
            return HttpResponseRedirect(
                redirect_to=f"/login?next={request.get_full_path()}"
            )
