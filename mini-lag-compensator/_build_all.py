import os
base = r"F:
ano-everything\mini-automation-theory. mini-classical-compensator\mini-lag-compensator"

def wf(rp, content):
    path = os.path.join(base, rp)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    lc = content.count(chr(10))
    print(f"OK {rp} ({lc} lines)")
