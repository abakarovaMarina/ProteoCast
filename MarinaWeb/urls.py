# MarinaWeb/urls.py
from django.contrib import admin
from django.urls import path, include
from browser.views import search_view, contact_us, documentation, citation, upload_file, drosophiladb, job_running

urlpatterns = [
    path('admin/', admin.site.urls),
    path('', search_view, name='home'),
    path('search/', include('browser.urls')),
    path('results/', include('browser.urls')),
    path('contact/', contact_us, name='contact_us'),
    path('documentation/', documentation, name='documentation'),
    path('citations/', citation, name='citation'),
    path('upload/', upload_file, name='upload'),
    path('job_running/<str:job_id>/', job_running, name='job_running'),
    path('drosophiladb/',drosophiladb, name='drosophilasb'),
]
