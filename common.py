from coffea.nanoevents import DelphesSchema

DelphesSchema.mixins.update({
    "GenParticle": "Particle",
    "GenCandidate": "Particle",
    "ParticleFlowCandidate": "Particle",
    "DarkHadronCandidate": "Particle",
    "FatJet": "Jet",
    "GenFatJet": "Jet",
    "DarkHadronJet": "Jet",
})

def load_sample(sample,helper=None,schema=DelphesSchema):
    from coffea.nanoevents import NanoEventsFactory
    path = f'models/{sample["model"]}'
    if helper is None:
        from svjHelper import svjHelper
        sample["helper"] = svjHelper.build(f'{path}/config.py')
    else:
        sample["helper"] = helper
    metadict = sample["helper"].metadata()
    metadict["dataset"] = sample["name"]
    sample["events"] = load_events(f'{path}/events.root',schema=schema,metadict=metadict)

def load_events(filename,schema=DelphesSchema,metadict=None):
    from coffea.nanoevents import NanoEventsFactory
    return NanoEventsFactory.from_root(
        file=filename,
        treepath="Delphes",
		schemaclass=schema,
		metadata=metadict,
	).events()
