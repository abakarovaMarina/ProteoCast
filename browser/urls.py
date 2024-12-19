from django.urls import path
from .views import search_view, results_view, download_folder, upload_file, drosophiladb, job_running

urlpatterns = [
    path('results/download/<str:fbpp_id>/', download_folder, name='download_file'),
    path('search/', search_view, name='search'),
    path('results/', results_view, name='results'),
    path('job_running/<str:job_id>/', job_running, name='job_running'),
    path('upload/', upload_file, name='upload_file'),
    path('drosophiladb/',drosophiladb, name='drosophiladb'),
]
