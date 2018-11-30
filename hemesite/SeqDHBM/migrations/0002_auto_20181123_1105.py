# Generated by Django 2.1.3 on 2018-11-23 10:05

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('SeqDHBM', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='job',
            name='full_hbm_analysis',
            field=models.TextField(blank=True),
        ),
        migrations.AlterField(
            model_name='sequence',
            name='status',
            field=models.CharField(choices=[('1', 'Received'), ('2', 'Processed'), ('3', 'Failed'), ('4', 'Skipped')], max_length=1),
        ),
        migrations.AlterField(
            model_name='sequence',
            name='submittedas',
            field=models.CharField(choices=[('M', 'Manual Input'), ('I', 'PDB Id'), ('F', 'Fasta File'), ('P', 'PDB File')], max_length=1),
        ),
    ]