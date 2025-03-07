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
    related_name = "compounds"
    # identifier = models.TextField(
    #     null=True, blank=True
    # )  # the identifier field is not used anywhere and by nooone. it is probaly some producer's identifier which is not
    # relevant
    name = models.TextField()
    structure = models.TextField(null=True, blank=True)
    # library = models.ForeignKey(
    #     CompoundLibrary,
    #     null=True,
    #     on_delete=models.SET_NULL,
    #     related_name=related_name,
    #     blank=True,
    # )
    data = models.JSONField(null=True)
    manufacturer = models.TextField(null=True, blank=True)
    item_number = models.TextField(null=True, blank=True)

    def __str__(self):
        return f"{self.name}"

    @property
    def structure_image(self):
        if not self.structure:
            return None
        mol = Chem.MolFromSmiles(self.structure)
        img = Draw.MolToImage(mol)
        buffer = BytesIO()
        img.save(buffer, format="PNG")
        return f"data:image/png;base64,{b64encode(buffer.getvalue()).decode()}"
