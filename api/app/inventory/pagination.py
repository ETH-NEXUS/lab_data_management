from rest_framework.pagination import PageNumberPagination


class InventoryStockPagination(PageNumberPagination):
    """
    Page-number pagination for inventory stock tables.

    Query examples:
    - `?page=1`
    - `?page=2&page_size=50`
    """

    page_size = 50
    page_size_query_param = "page_size"
    max_page_size = 200


class InventoryHistoryPagination(PageNumberPagination):
    """
    Page-number pagination for inventory history records.

    Query examples:
    - `?page=1&page_size=5`
    - `?page=2&page_size=50`
    """

    page_size = 50
    page_size_query_param = "page_size"
    max_page_size = 200
