# MarinaWeb/urls.py
from django.contrib import admin
from django.urls import path, include
from browser.views import search_view 

urlpatterns = [
    path('admin/', admin.site.urls),
    path('', search_view, name='home'),  
    path('search/', include('browser.urls')),  
    path('results/', include('browser.urls')), 
]
