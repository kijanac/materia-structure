import qcio

p = qcio.Packmol()

settings = qcio.Settings()
settings.run.executable = "/home/kijana/sw/packmol-20.2.2/packmol"
settings.run.inp = "pack.in"
settings.run.log = "pack.out"
settings.run.out = "pack.xyz"
settings.packmol.tolerance = 2.0
settings.packmol.structures = ["h2o.xyz", "h2o.xyz"]
settings.packmol.counts = [30, 100]
settings.packmol.constraints = [
    "inside box 0. 0. 0. 10. 10. 10.",
    "outside box 0. 0. 0. 10. 10. 10.",
]
settings.packmol.nloop = 200

p(settings, 1)
