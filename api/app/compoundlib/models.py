from django.db import models
from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO
from base64 import b64encode


class CompoundLibrary(models.Model):
    name = models.CharField(max_length=50, unique=True, verbose_name="compound library")
    file_name = models.CharField(max_length=255, null=True, blank=True)

    def __str__(self):
        return f"{self.name}"

    class Meta:
        verbose_name_plural = "compound libraries"
        ordering = ("name",)


class Compound(models.Model):
    identifier = models.TextField(unique=True)
    name = models.TextField()
    structure = models.TextField()
    library = models.ForeignKey(CompoundLibrary, null=True, on_delete=models.RESTRICT)
    data = models.JSONField(null=True)

    def __str__(self):
        return f"{self.name}"

    @property
    def structure_image(self):
        mol = Chem.MolFromSmiles(self.structure)
        img = Draw.MolToImage(mol)
        buffer = BytesIO()
        img.save(buffer, format="PNG")
        return f"data:image/png;base64,{b64encode(buffer.getvalue()).decode()}"
