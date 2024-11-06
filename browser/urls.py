from django.urls import path
from .views import search_view, results, download_file

urlpatterns = [
    path('', search_view, name='search'),
    path('search/', search_view, name='search'),
    path('results/', results, name='results'),
    path('download/<str:fbpp_id>/', download_file, name='download_file'),
]