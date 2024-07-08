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
