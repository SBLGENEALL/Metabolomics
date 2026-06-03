"""
setup_env.py — 초기 환경 설정 (처음 한 번만)
사용: python setup_env.py
"""
import subprocess, sys, os, zipfile, io

ROOT = os.path.dirname(os.path.abspath(__file__))
print("="*55); print("  CHO_METABOLOMICS 환경 설정"); print("="*55)

PKGS = ["cobra","pandas","numpy","scipy","matplotlib",
        "seaborn","openpyxl","requests","scikit-learn","escher"]
missing = []
for pkg in PKGS:
    imp = {"scikit-learn":"sklearn","escher":"escher"}.get(pkg,pkg)
    try: __import__(imp); print(f"  ✔ {pkg}")
    except ImportError: missing.append(pkg); print(f"  ✘ {pkg}")
if missing:
    print(f"\n  설치 중: {missing}")
    subprocess.check_call([sys.executable,"-m","pip","install"]+missing)

MODEL_DIR = os.path.join(ROOT,"model","iCHO3K")
os.makedirs(MODEL_DIR, exist_ok=True)
prod = next((f for f in os.listdir(MODEL_DIR) if "prod" in f and f.endswith(".json")),None)
if prod:
    print(f"\n  ✔ 모델: {prod}")
else:
    print("\n  iCHO3K 다운로드 중 (LewisLabUCSD/iCHO3K)...")
    try:
        import requests
        r = requests.get("https://github.com/LewisLabUCSD/iCHO3K/archive/refs/heads/main.zip",timeout=180)
        r.raise_for_status()
        with zipfile.ZipFile(io.BytesIO(r.content)) as z:
            for m in z.namelist():
                if "iCHO3K/Model/" in m and m.endswith(".json"):
                    fname = os.path.basename(m)
                    with z.open(m) as s, open(os.path.join(MODEL_DIR,fname),"wb") as d:
                        d.write(s.read())
                    print(f"  ✔ {fname}")
    except Exception as e:
        print(f"  ✘ {e}")
        print("  수동: github.com/LewisLabUCSD/iCHO3K → JSON → model/iCHO3K/")

sys.path.insert(0,ROOT)
from src.config import MODEL_PATH,SOWA_MMC1,AA20_EXCEL
print("\n  [파일 확인]")
for k,v in [("모델",MODEL_PATH),("Sowa",SOWA_MMC1),("20AA",AA20_EXCEL)]:
    print(f"  {'✔' if v and os.path.exists(v) else '✘'} {k}: {v or '없음'}")
print("\n  완료! python run_pipeline.py --dataset practice_20aa --rate_days 7,10")
