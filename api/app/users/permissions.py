from rest_framework import permissions


class IsOwnUser(permissions.BasePermission):
    """
    Object-level permission to only allow to get the user object of himself.
    """

    def has_object_permission(self, request, view, obj):
        return obj.id == request.user.id
