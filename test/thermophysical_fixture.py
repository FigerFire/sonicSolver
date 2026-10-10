"""Complete old regression declarations in disposable test copies only.

These fixtures were frozen with gamma=1.4/R=287.05 and numerics.transport
coefficients. This helper makes those specific historical choices explicit;
it is never imported by the solver or applied to production case directories.
"""
from pathlib import Path
import re

def declare_historical_thermophysical_fixture(case):
    case=Path(case)
    path=case/'models/thermoDynamics.yaml'
    if not path.exists():
        return
    text=path.read_text()
    numerics=(case/'solvers/numerics.yaml').read_text()
    def coefficient(name, historical):
        match=re.search(r'^\s+'+name+r':\s*([^\n#]+)',numerics,re.M)
        return match.group(1).strip() if match else historical
    mu=coefficient('mu','0.0');Pr=coefficient('Pr','0.72')
    if 'properties: {}' in text:
        if 'equationOfState: perfectGas' not in text:
            raise AssertionError('Test fixture is not historical PerfectGas: '+str(path))
        text=text.replace('properties: {}','properties:\n  equationOfState:\n    gamma: 1.4\n    R: 287.05\n  transport:\n    mu: '+mu+'\n    Pr: '+Pr)
    elif 'equationOfState: rhoConst' in text and not re.search(r'^  transport:[ \t]*$',text,re.M):
        text+='\n  transport:\n    mu: '+mu+'\n    Pr: '+Pr+'\n'
    path.write_text(text)
