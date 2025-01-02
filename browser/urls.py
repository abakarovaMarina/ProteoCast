from django.urls import path
from .views import search_view, results_view, download_folder, upload_file, drosophiladb, job_running, serve_file, check_job_status, results_job

urlpatterns = [
    path('search/', search_view, name='search'),
    path('results/', results_view, name='results'),
    path('results_job/<str:job_id>/', results_job, name='results_job'),
    path('drosophiladb/',drosophiladb, name='drosophiladb'),
    path('data/<str:folder>/<str:filename>/', serve_file, name='serve_file'),
    path('download_folder/<str:fbpp_id>/', download_folder, name='download_folder'),
    path('upload_file/', upload_file, name='upload_file'),
    path('job_running/<str:job_id>/', job_running, name='job_running'),
    path('check_job_status/', check_job_status, name='check_job_status'), 

]

