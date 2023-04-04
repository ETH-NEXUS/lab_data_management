# Generated by Django 4.1.7 on 2023-04-03 11:03

import core.basemodels
import core.models
import django.contrib.postgres.fields
import django.core.validators
from django.db import migrations, models
import django.db.models.deletion
import django.utils.timezone


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        ('platetemplate', '0004_alter_platetemplate_category'),
        ('compoundlib', '0008_alter_compoundlibrary_options'),
    ]

    operations = [
        migrations.CreateModel(
            name='PlateDetail',
            fields=[
                ('id', models.BigIntegerField(primary_key=True, serialize=False)),
                ('num_wells', models.IntegerField()),
                ('measurement_labels', django.contrib.postgres.fields.ArrayField(base_field=models.TextField(blank=True, null=True), size=None)),
                ('stats', core.models.DictField()),
                ('overall_stats', core.models.DictField()),
            ],
            options={
                'db_table': 'core_platedetail',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='WellDetail',
            fields=[
                ('id', models.BigIntegerField(primary_key=True, serialize=False)),
                ('plate_id', models.BigIntegerField()),
                ('hr_position', models.CharField(max_length=10)),
                ('initial_amount', models.FloatField(blank=True, null=True)),
                ('withdrawal', models.FloatField(blank=True, null=True)),
                ('amount', models.FloatField(blank=True, null=True)),
                ('compounds', django.contrib.postgres.fields.ArrayField(base_field=models.TextField(blank=True, null=True), size=None)),
                ('measurements', core.models.DictField()),
            ],
            options={
                'db_table': 'core_welldetail',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='Experiment',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created_at', models.DateTimeField(default=django.utils.timezone.now, editable=False)),
                ('modified_at', core.basemodels.AutoDateTimeField(editable=False)),
                ('name', models.CharField(max_length=50)),
                ('description', models.TextField(blank=True, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='Location',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created_at', models.DateTimeField(default=django.utils.timezone.now, editable=False)),
                ('modified_at', core.basemodels.AutoDateTimeField(editable=False)),
                ('name', models.CharField(max_length=50, verbose_name='location')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='MeasurementFeature',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('abbrev', models.CharField(max_length=20, unique=True)),
                ('name', models.CharField(blank=True, max_length=50, null=True, verbose_name='measurement')),
                ('unit', models.CharField(blank=True, max_length=10, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='Plate',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created_at', models.DateTimeField(default=django.utils.timezone.now, editable=False)),
                ('modified_at', core.basemodels.AutoDateTimeField(editable=False)),
                ('barcode', models.CharField(db_index=True, max_length=50, unique=True, validators=[django.core.validators.RegexValidator('[^\\s]+', 'Plate barcode must not contain strings!')])),
            ],
            options={
                'ordering': ('-id',),
            },
        ),
        migrations.CreateModel(
            name='PlateDimension',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=50, verbose_name='plate dimension')),
                ('rows', models.PositiveIntegerField(validators=[django.core.validators.MinValueValidator(1)])),
                ('cols', models.PositiveIntegerField(validators=[django.core.validators.MinValueValidator(1)])),
            ],
        ),
        migrations.CreateModel(
            name='Project',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created_at', models.DateTimeField(default=django.utils.timezone.now, editable=False)),
                ('modified_at', core.basemodels.AutoDateTimeField(editable=False)),
                ('name', models.CharField(max_length=50, unique=True)),
                ('description', models.TextField(blank=True, null=True)),
                ('harvest_id', models.IntegerField(blank=True, null=True)),
                ('harvest_notes', models.TextField(blank=True, null=True)),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='Sample',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created_at', models.DateTimeField(default=django.utils.timezone.now, editable=False)),
                ('modified_at', core.basemodels.AutoDateTimeField(editable=False)),
                ('name', models.CharField(max_length=50, verbose_name='sample')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='Well',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created_at', models.DateTimeField(default=django.utils.timezone.now, editable=False)),
                ('modified_at', core.basemodels.AutoDateTimeField(editable=False)),
                ('position', models.PositiveIntegerField(db_index=True)),
                ('status', models.TextField(blank=True, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='WellType',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(db_index=True, max_length=50)),
                ('description', models.TextField()),
            ],
        ),
        migrations.CreateModel(
            name='WellWithdrawal',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created_at', models.DateTimeField(default=django.utils.timezone.now, editable=False)),
                ('modified_at', core.basemodels.AutoDateTimeField(editable=False)),
                ('amount', models.FloatField()),
                ('target_well', models.ForeignKey(null=True, on_delete=django.db.models.deletion.SET_NULL, related_name='donors', to='core.well')),
                ('well', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='withdrawals', to='core.well')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='WellCompound',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('amount', models.FloatField(default=0, validators=[django.core.validators.MinValueValidator(0)])),
                ('compound', models.ForeignKey(on_delete=django.db.models.deletion.RESTRICT, related_name='well_compounds', to='compoundlib.compound')),
                ('well', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='well_compounds', to='core.well')),
            ],
            options={
                'unique_together': {('well', 'compound')},
            },
        ),
        migrations.AddField(
            model_name='well',
            name='compounds',
            field=models.ManyToManyField(through='core.WellCompound', to='compoundlib.compound'),
        ),
        migrations.AddField(
            model_name='well',
            name='plate',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='wells', to='core.plate'),
        ),
        migrations.AddField(
            model_name='well',
            name='sample',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.RESTRICT, related_name='wells', to='core.sample'),
        ),
        migrations.AddField(
            model_name='well',
            name='type',
            field=models.ForeignKey(default=1, on_delete=django.db.models.deletion.RESTRICT, to='core.welltype'),
        ),
        migrations.CreateModel(
            name='PlateMapping',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created_at', models.DateTimeField(default=django.utils.timezone.now, editable=False)),
                ('modified_at', core.basemodels.AutoDateTimeField(editable=False)),
                ('mapping_file', models.FileField(null=True, upload_to='')),
                ('from_column', models.CharField(max_length=50, null=True)),
                ('to_column', models.CharField(max_length=50, null=True)),
                ('amount_column', models.CharField(max_length=50, null=True)),
                ('delimiter', models.CharField(default=',', max_length=1, null=True)),
                ('quotechar', models.CharField(default='"', max_length=1, null=True)),
                ('amount', models.FloatField(default=None, null=True, validators=[django.core.validators.MinValueValidator(0)])),
                ('evaluation', models.TextField(blank=True, null=True)),
                ('source_plate', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='mapped_to_plates', to='core.plate')),
                ('target_plate', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='mapped_from_plates', to='core.plate')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.AddField(
            model_name='plate',
            name='dimension',
            field=models.ForeignKey(default=None, null=True, on_delete=django.db.models.deletion.RESTRICT, to='core.platedimension'),
        ),
        migrations.AddField(
            model_name='plate',
            name='experiment',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='plates', to='core.experiment'),
        ),
        migrations.AddField(
            model_name='plate',
            name='library',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='plates', to='compoundlib.compoundlibrary'),
        ),
        migrations.AddField(
            model_name='plate',
            name='template',
            field=models.OneToOneField(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='plate', to='platetemplate.platetemplate'),
        ),
        migrations.CreateModel(
            name='MeasurementAssignment',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created_at', models.DateTimeField(default=django.utils.timezone.now, editable=False)),
                ('modified_at', core.basemodels.AutoDateTimeField(editable=False)),
                ('status', models.CharField(default='pending', max_length=50)),
                ('filename', models.TextField()),
                ('measurement_file', models.FileField(null=True, upload_to='')),
                ('plate', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='assignments', to='core.plate')),
            ],
            options={
                'abstract': False,
            },
        ),
        migrations.CreateModel(
            name='Measurement',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created_at', models.DateTimeField(default=django.utils.timezone.now, editable=False)),
                ('modified_at', core.basemodels.AutoDateTimeField(editable=False)),
                ('value', models.FloatField()),
                ('label', models.CharField(blank=True, max_length=50, null=True)),
                ('identifier', models.CharField(blank=True, max_length=20, null=True)),
                ('measured_at', models.DateTimeField(blank=True, null=True)),
                ('feature', models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.RESTRICT, related_name='measurements', to='core.measurementfeature')),
                ('measurement_assignment', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='measurements', to='core.measurementassignment')),
                ('well', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='measurements', to='core.well')),
            ],
        ),
        migrations.AddField(
            model_name='experiment',
            name='project',
            field=models.ForeignKey(on_delete=django.db.models.deletion.RESTRICT, related_name='experiments', to='core.project'),
        ),
        migrations.CreateModel(
            name='BarcodeSpecification',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('created_at', models.DateTimeField(default=django.utils.timezone.now, editable=False)),
                ('modified_at', core.basemodels.AutoDateTimeField(editable=False)),
                ('prefix', models.CharField(max_length=100)),
                ('number_of_plates', models.IntegerField(blank=True, null=True)),
                ('sides', django.contrib.postgres.fields.ArrayField(base_field=models.CharField(max_length=20), blank=True, null=True, size=None)),
                ('experiment', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='barcode_specifications', to='core.experiment')),
            ],
            options={
                'ordering': ['id'],
            },
        ),
        migrations.AlterUniqueTogether(
            name='well',
            unique_together={('plate', 'position')},
        ),
        migrations.AddConstraint(
            model_name='plate',
            constraint=models.CheckConstraint(check=models.Q(models.Q(('experiment__isnull', True), ('library__isnull', True), ('template__isnull', True)), models.Q(('experiment__isnull', False), ('library__isnull', True), ('template__isnull', True)), models.Q(('experiment__isnull', True), ('library__isnull', False), ('template__isnull', True)), models.Q(('experiment__isnull', True), ('library__isnull', True), ('template__isnull', False)), _connector='OR'), name='check_only_library_or_experiment_or_template'),
        ),
        migrations.AddConstraint(
            model_name='plate',
            constraint=models.UniqueConstraint(fields=('barcode', 'experiment'), name='unique_barcode_experiment'),
        ),
        migrations.AlterUniqueTogether(
            name='measurement',
            unique_together={('well', 'label', 'measured_at')},
        ),
        migrations.AlterUniqueTogether(
            name='experiment',
            unique_together={('name', 'project')},
        ),
    ]
