# Generated by Django 2.1.3 on 2018-11-23 16:33

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('SeqDHBM', '0006_auto_20181123_1346'),
    ]

    operations = [
        migrations.AddField(
            model_name='sequence',
            name='warnings_hbm',
            field=models.TextField(blank=True),
        ),
    ]
