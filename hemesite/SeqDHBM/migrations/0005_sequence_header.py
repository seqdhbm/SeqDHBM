# Generated by Django 2.1.3 on 2018-11-23 12:22

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('SeqDHBM', '0004_auto_20181123_1216'),
    ]

    operations = [
        migrations.AddField(
            model_name='sequence',
            name='header',
            field=models.CharField(default='Forgot to save the sequence header', max_length=300),
            preserve_default=False,
        ),
    ]
