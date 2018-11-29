from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('runworkflow/', views.hemewf, name='hemewf'),
    path('<int:job_id>/', views.show_result, name='result'),
    # path(<route>, <view>, <kwargs>, <name>)
    # <route> contains a URL pattern. Django starts at the first pattern in urlpatterns
    #  and makes its way down the list, comparing the requested URL against each pattern until it finds one that matches
    # <view> When Django finds a matching pattern, it calls the specified view function with an HttpRequest object
    # as the first argument and any “captured” values from the route as keyword arguments.
    # <kwargs> parameters
    # <name> Naming your URL lets you refer to it unambiguously from elsewhere in Django, especially from within templates
]
