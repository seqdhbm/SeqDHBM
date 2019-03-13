from django.urls import path
from django.contrib.staticfiles.urls import staticfiles_urlpatterns

from . import views

# requests are usually handled by the module hemesite/SeqDHBM/views.py
urlpatterns = [
    # handles requests to the main page
    path('', views.index, name='index'),
    # Handles requests for results without the password (hide information)
    path('<int:job_id>/', views.show_result, name='result'),
    # Handles requests for results
    path('<int:job_id>/<passw>/', views.show_result, name='result'),
    # Shows the full analysis
    path('analysis/<int:job_id>/<passw>/', views.show_analysis, name='analysis'),

    path('runworkflow/<int:job_id>/<passw>/', views.hemewf, name='hemewf'),

    path('analysis/<int:job_id>/', views.show_analysis, name='analysis'),
    # path(<route>, <view>, <kwargs>, <name>)
    # <route> contains a URL pattern. Django starts at the first pattern in urlpatterns
    #  and makes its way down the list, comparing the requested URL against each pattern until it finds one that matches
    # <view> When Django finds a matching pattern, it calls the specified view function with an HttpRequest object
    # as the first argument and any “captured” values from the route as keyword arguments.
    # <kwargs> parameters
    # <name> Naming your URL lets you refer to it unambiguously from elsewhere in Django, especially from within templates
]

urlpatterns += staticfiles_urlpatterns()