# coding: utf-8

"""
Django settings for hemesite project.

For more information on this file, see
https://docs.djangoproject.com/en/2.1/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/2.1/ref/settings/
"""

import os

from . import local_settings as cfg

# It's secret o.O
SECRET_KEY = cfg.SECRET_KEY

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Celery settings

CELERY_BROKER_URL = 'amqp://guest:guest@localhost:5672//'
#: Only add pickle to this list if your broker is secured
#: from unwanted access (see userguide/security.html)
CELERY_ACCEPT_CONTENT = ['json']
CELERY_RESULT_BACKEND = 'db+sqlite:///results.sqlite'
CELERY_TASK_SERIALIZER = 'json'
CELERY_TIMEZONE = "UTC"
CELERY_ENABLE_UTC = True
CELERY_IMPORTS = ['seqdhbm']


# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/2.1/howto/deployment/checklist/

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = cfg.DEBUG


ALLOWED_HOSTS = ["131.220.127.75", "131.220.127.81", "127.0.0.1", "localhost", "rms1705102"]

# This will be used to send links to the user access the results
SITE_DOMAIN = cfg.SITE_DOMAIN

# Application definition

INSTALLED_APPS = [
    'SeqDHBM.apps.SeqdhbmConfig',
    'django.contrib.admin',
    'django_celery_beat',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
]

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

ROOT_URLCONF = 'hemesite.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

WSGI_APPLICATION = 'hemesite.wsgi.application'

# Database
# https://docs.djangoproject.com/en/2.1/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME': os.path.join(BASE_DIR, 'db.sqlite3'),
    }
}

# Password validation
# https://docs.djangoproject.com/en/2.1/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]


# Internationalization
# https://docs.djangoproject.com/en/2.1/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'Europe/Berlin'

USE_I18N = True

USE_L10N = True

USE_TZ = True

# Emails
# Mail is sent using the SMTP host and port specified in the
EMAIL_HOST = "smtp.gmail.com"
EMAIL_PORT = 465  # 587 for TLS
EMAIL_USE_TLS = False
EMAIL_USE_SSL = True
# control whether a secure connection is used

EMAIL_HOST_USER = "seqdhbm@gmail.com"
# SECURITY WARNING: keep the secret key used in production secret!
EMAIL_HOST_PASSWORD = cfg.EMAIL_HOST_PASSWORD

# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/2.1/howto/static-files/

STATIC_URL = '/static/'

STATICFILES_DIRS = [
    os.path.join(BASE_DIR, "static"),
]
# Security measures

# Change SECURE_PROXY_SSL_HEADER if the server is behind a proxy. Check the
# link for more information:
# https://docs.djangoproject.com/en/2.1/ref/settings/#std:setting-SECURE_PROXY_SSL_HEADER
# SECURE_PROXY_SSL_HEADER = None

SECURE_SSL_REDIRECT = False
# SECURE_HSTS_SECONDS = 20

MAX_FILE_SIZE = 2097152  # 2MB - There are not many fasta files that large

# DEPLOYMENT
WSGI_APPLICATION
