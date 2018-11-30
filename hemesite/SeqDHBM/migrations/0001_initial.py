# Generated by Django 2.1.3 on 2018-11-22 11:23

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Job',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('submittedby', models.EmailField(blank=True, max_length=254)),
                ('submission_date', models.DateTimeField(auto_now_add=True)),
            ],
        ),
        migrations.CreateModel(
            name='Result',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('coord_atom', models.CharField(max_length=6)),
                ('ninemer', models.CharField(max_length=9)),
                ('net_charge', models.SmallIntegerField()),
                ('disulfide_possible', models.BooleanField()),
            ],
        ),
        migrations.CreateModel(
            name='Sequence',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('seqchain', models.TextField()),
                ('submittedas', models.CharField(choices=[('M', 'Manual Input'), ('F', 'Fasta File'), ('P', 'PDB File'), ('I', 'PDB Id')], max_length=1)),
                ('mode', models.CharField(choices=[('W', 'wesa'), ('S', 'structure')], max_length=1)),
                ('status', models.CharField(choices=[('1', 'Received'), ('2', 'Processed'), ('3', 'Failed')], max_length=1)),
                ('jobnum', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='SeqDHBM.Job', verbose_name='Part of the job')),
            ],
        ),
        migrations.AddField(
            model_name='result',
            name='sequence',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='SeqDHBM.Sequence', verbose_name='Belongs to a sequence'),
        ),
    ]