
class LorentzFactor(object):
    def transform(self, prim_selection):
        v1 = prim_selection['v1']
        v2 = prim_selection['v2']
        v3 = prim_selection['v3']
        return (1.0 - (v1**2 + v2**2 + v3**2))**-0.5

