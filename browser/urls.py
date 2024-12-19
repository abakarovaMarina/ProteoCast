from django.urls import path
from .views import search_view, results_view, download_folder, upload_file, drosophiladb, job_running, serve_file

urlpatterns = [
    path('search/', search_view, name='search'),
    path('results/', results_view, name='results'),
    path('job_running/', job_running, name='job_running'),
    path('upload/', upload_file, name='upload_file'),
    path('drosophiladb/',drosophiladb, name='drosophiladb'),
    path('data/<str:folder>/<str:filename>/', serve_file, name='serve_file'),
    path('download_folder/<str:fbpp_id>/', download_folder, name='download_folder'),
]
