"""Korean-language versions of all five figures.

Each figure script is read, its display strings replaced, the font switched to a
CJK face, and the output written as Figure{n}_ko.png/pdf. The computation is
untouched: only text is translated.
"""
import re, pathlib, subprocess, sys

SRC = pathlib.Path(__file__).resolve().parent
TMP = pathlib.Path("/tmp/claude-0/-home-claude/5bca9a6a-66b9-5ed3-a5ad-29e20b71c1e0/scratchpad/kofig")
TMP.mkdir(parents=True, exist_ok=True)

FONT_HEADER = '''
import matplotlib, matplotlib.font_manager as _fm
_ko = sorted({f.name for f in _fm.fontManager.ttflist if "Noto Sans CJK" in f.name})
matplotlib.rcParams["font.family"] = _ko[0]
matplotlib.rcParams["axes.unicode_minus"] = False
'''

COMMON = {
 "not observable": "관측 불가",
 "not specific": "비특이적",
 "dataset-specific": "데이터셋 한정",
 "generalizable": "일반화 가능",
 "Ward": "Ward",
 "k-means": "k-means",
}

REPL = {
"fig1": {
 "A reported boundary enters the audit": "보고된 경계가 감사에 들어온다",
 "reported\\nboundary": "보고된\\n경계",
 "G1  observability": "G1  관측가능성",
 "G2  specificity": "G2  특이성",
 "G3  independence": "G3  독립성",
 "G4  donor stability": "G4  공여자 안정성",
 "G5  region / cohort": "G5  영역 · 코호트",
 "which markers clear the detection\\nfloor and the ambient band?":
   "어느 마커가 검출 하한과\\n앰비언트 대역을 넘는가?",
 "is it absent from matched random\\npanels, scored on the same statistic?":
   "같은 통계량으로 채점했을 때\\n맞춘 무작위 패널에 없는가?",
 "do the agreeing components carry\\nindependent evidence?":
   "동의하는 성분들이 서로\\n별개의 증거인가?",
 "does it survive leave-one-donor-out\\nand a donor bootstrap?":
   "공여자 1인 제외와 부트스트랩을\\n견디는가?",
 "is it recovered where composition\\nallows the same question?":
   "세포 구성이 같은 질문을 허용하는\\n곳에서 재현되는가?",
 "not specific\\n— a property of the grid": "비특이적\\n— 격자의 성질",
 "unevaluable\\n— agreement is redundancy": "평가 불가\\n— 일치는 중복",
 "unevaluable\\n— carried by one donor": "평가 불가\\n— 한 사람이 떠받침",
 "the earliest failing gate sets the verdict; every gate outcome is reported":
   "가장 먼저 실패하는 관문이 판정을 정하며, 모든 관문의 결과를 보고한다",
 "The three claims audited here, and the gate at which each stops":
   "여기서 감사한 세 주장과 각각이 멈추는 관문",
 "CPS 0.4–0.6\\n(atlas authors)": "CPS 0.4–0.6\\n(아틀라스 저자)",
 "CPS 0.207\\n(this panel)": "CPS 0.207\\n(이 패널)",
 "HMOX1\\n(largest effect)": "HMOX1\\n(최대 효과)",
 "G1 pass\\nG2  p = 0.079": "G1 통과\\nG2  p = 0.079",
 "G1 pass on restriction\\nG2  p = 0.753": "G1 제한 후 통과\\nG2  p = 0.753",
 "G1 fail\\n25 of 84 donors": "G1 실패\\n84명 중 25명",
 "not specific\\n(marginal)": "비특이적\\n(경계선상)",
 "No claim examined here reaches \\u201cgeneralizable\\u201d. That the audit certifies ":
   "여기서 검토한 어떤 주장도 「일반화 가능」에 이르지 못한다. 이 감사가 합성 패널 16개 중 ",
 "planted boundaries in 16 of 16 synthetic panels (Figure 2) is what\\nseparates this result ":
   "16개에서 심어 둔 경계를 인증한다는 사실(그림 2)이\\n이 결과를 모든 것을 기각하는 ",
 "from a procedure that rejects everything. A null at G1 is not evidence of no change — only ":
   "절차와 구별한다. G1에서의 음성은 변화가 없다는 증거가 아니라 ",
 "that the assay cannot see it.": "분석이 그것을 볼 수 없다는 뜻일 뿐이다.",
},
"fig2": {
 "real boundary": "실제 경계", "global offset": "전역 이동", "null panel": "귀무 패널",
 "monotone ramp": "단조 경사", "single-marker decoy": "단일 마커 미끼",
 "single-donor artefact": "단일 공여자 인공물",
 "published gates (3)": "발표된 관문 (3)", "+ leave-one-donor-out": "+ 공여자 1인 제외",
 "panels certified (of {NREP})": "인증된 패널 수 (총 {NREP})",
 "what the audit certifies": "감사가 인증하는 것",
 "% of panels passing each gate": "관문별 통과 패널 비율 (%)",
 "speci-\\nficity": "특이성", "marker\\nLOO": "마커\\n제외",
 "donor\\nLOO": "공여자\\n제외",
 "real\\nboundary": "실제\\n경계", "single-donor\\nartefact": "단일 공여자\\n인공물",
 "largest shift when one\\ndonor is removed (CPS)": "공여자 1인 제외 시\\n최대 경계 이동 (CPS)",
 "gate 0.05": "관문 0.05",
 "what leave-one-donor-out sees": "공여자 1인 제외가 보는 것",
},
"fig3": {
 "equal-width bins are\\nnot equal-sized": "같은 너비의 bin이\\n같은 크기는 아니다",
 "CPS bin centre": "CPS bin 중심", "donors in bin": "bin 안의 공여자 수",
 "boundary": "경계",
 "detection rate in astrocytes (log scale)": "성상교세포 검출률 (로그 척도)",
 "floor 0.02": "하한 0.02", "ambient": "앰비언트",
 "below floor": "하한 미만", "too few donors": "공여자 부족",
 "at or below ambient": "앰비언트 이하", "observable": "관측 가능",
 "cannot be adjudicated": "는 판정 불가",
 "dropping the 23 changes\\nnothing": "23개를 빼도\\n달라지지 않는다",
 "Ward first boundary (CPS)": "Ward 1차 경계 (CPS)",
 "all {int(res.n_panel)}\\nmarkers": "마커\\n{int(res.n_panel)}개 전부",
 "{int(res.n_observable)} observable\\nonly": "관측 가능\\n{int(res.n_observable)}개만",
 "boundary with one donor removed (CPS)": "공여자 1인 제외 시 경계 (CPS)",
 "{int((~inside).sum())} of {len(L)} donors move it":
   "{len(L)}명 중 {int((~inside).sum())}명이 경계를 움직인다",
 "donor bootstrap (CPS)": "공여자 부트스트랩 (CPS)", "resamples": "재표본 수",
 "but the interval is wide": "다만 구간은 넓다",
 "median {np.median(bv):.3f}\\n95% CI {lo:.3f}–{hi:.3f}\\n":
   "중앙값 {np.median(bv):.3f}\\n95% CI {lo:.3f}–{hi:.3f}\\n",
 "{np.mean(np.abs(bv-BP)<=TOL)*100:.0f}% within 0.05":
   "{np.mean(np.abs(bv-BP)<=TOL)*100:.0f}%가 0.05 이내",
 "{int(res.n_panel - res.n_observable)} of {int(res.n_panel)} panel markers ":
   "패널 마커 {int(res.n_panel)}개 중 {int(res.n_panel - res.n_observable)}개",
},
"fig4": {
 "PC1 breakpoint  (the null as published)": "PC1 변화점  (발표된 귀무분포)",
 "Ward boundary  (like for like)": "Ward 경계  (동일 통계량)",
 "observed\\n0.207": "관측값\\n0.207",
 "breakpoint in 1,000 expression-matched random panels (CPS)":
   "발현량 일치 무작위 패널 1,000개의 변화점 (CPS)",
 "panels": "패널 수",
 "the same boundary, two nulls": "같은 경계, 두 개의 귀무분포",
 "Ward 0.207\\nvs PC1 null": "Ward 0.207\\nvs PC1 귀무",
 "Ward 0.207\\nvs Ward null": "Ward 0.207\\nvs Ward 귀무",
 "PC1 0.543\\nvs PC1 null": "PC1 0.543\\nvs PC1 귀무",
 "as published": "발표된 방식", "like for like": "동일 통계량",
 "fraction of random panels within 0.05 CPS": "0.05 CPS 이내 무작위 패널 비율",
 "the p-value is a property of the null": "p값은 귀무분포의 성질이다",
 "any random 3-labelling": "임의의 무작위 3-라벨링",
 "chance (1 of 8 positions)": "우연 (8개 위치 중 1개)",
 "SEA-AD grid\\n(100%)": "SEA-AD 격자\\n(100%)",
 "donors in the first bin (nine-bin grid, pure noise)": "첫 bin의 공여자 수 (9개 bin, 순수 잡음)",
 "first cut lands at the\\nfirst midpoint (%)": "첫 절단이 첫 중점에\\n떨어지는 비율 (%)",
 "on noise, the grid decides where the cut falls": "잡음에서는 격자가 절단 위치를 정한다",
 "agreement is redundancy, not support": "일치는 지지가 아니라 중복이다",
 "random panels returning the identical value (%)": "동일한 값을 돌려준 무작위 패널 (%)",
 "and {len(moved)} of 16 components move more than 0.05 CPS when the\\n":
   "그리고 같은 행렬을 1% 다른 척도로 다시 세우면 16개 성분 중\\n",
 "same matrix is rebuilt with a 1% change in scaling — including\\n":
   "{len(moved)}개가 0.05 CPS 넘게 움직인다 — 그중 하나가 분산 점프\\n",
 "variance jump (0.207 \\u2192 0.886), one of the four supporting the boundary":
   "(0.207 → 0.886)로, 경계를 지지하던 넷 중 하나다",
},
"fig5": {
 "reported 0.4–0.6": "보고된 0.4–0.6", "panel boundary": "패널 경계",
 "staging, scaled 0–1": "병기, 0–1 척도", "CPS": "CPS",
 "neuropathological staging": "신경병리 병기",
 "Ward null\\n(random panels)": "Ward 귀무분포\\n(무작위 패널)",
 "panel": "패널",
 "{nin} of 16": "{nin}/16", "{nour} of 16": "{nour}/16",   # 「16개 중 n개」는 영문보다 길어 두 라벨이 겹친다
 "breakpoint on the four staging axes (CPS)": "네 병기 축에서의 변화점 (CPS)",
 "staging escapes the forced position": "병기는 강제된 위치에서 벗어난다",
 "Thal (Aβ)": "Thal (Aβ)", "Braak (tau)": "Braak (타우)", "dementia": "치매",
 "CPS is not tau-privileged": "CPS는 타우에 특권적이지 않다",
 "Spearman $\\\\rho$ with CPS": "CPS와의 Spearman \u03c1",   # 수식($..$)을 쓰면 matplotlib이 문자열 전체를 수식 글꼴로 그려 한글이 깨진다
},
}

if __name__ == "__main__":
    for name in ["fig1", "fig2", "fig3", "fig4", "fig5"]:
        src = (SRC / f"{name}.py").read_text()
        # 임시 디렉터리에서 실행하므로 paths 모듈이 있는 원본 디렉터리를 먼저 올린다
        src = src.replace('sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))',
                          f'sys.path.insert(0, {str(SRC)!r})')
        src = src.replace('from paths import DATA, INPUT, FIG',
                          'from paths import DATA, INPUT, FIG' + FONT_HEADER)
        if name == "fig1":
            import re as _re
            src = _re.sub(r"No claim examined here reaches .{0,20}generalizable.{0,20}\. That the audit certifies ",
                          "여기서 검토한 어떤 주장도 「일반화 가능」에 이르지 못한다. 이 감사가 합성 패널 16개 중 ", src)
        if name == "fig4":
            import re as _re
            src = _re.sub(r"variance jump \(0\.207 .{0,12} 0\.886\), one of the four supporting the boundary",
                          "(0.207 → 0.886)로, 경계를 지지하던 넷 중 하나다", src)
        if name == "fig2":
            src = src.replace('"support", "speci-\\nficity"', '"지지도", "특이성"')
        table = dict(COMMON); table.update(REPL[name])
        hits = 0
        for a, b in sorted(table.items(), key=lambda t: -len(t[0])):
            if a.startswith("GL = ") or a.startswith("RAW:"):
                k = a[4:] if a.startswith("RAW:") else a
                if k in src: src = src.replace(k, b); hits += 1
                continue
            needle = f'"{a}"'
            if needle in src:
                src = src.replace(needle, f'"{b}"'); hits += 1
        n = name[-1]
        src = src.replace(f'FIG + "Figure{n}.png"', f'FIG + "Figure{n}_ko.png"')
        src = src.replace(f'FIG + "Figure{n}.pdf"', f'FIG + "Figure{n}_ko.pdf"')
        out = TMP / f"{name}_ko.py"; out.write_text(src)
        r = subprocess.run([sys.executable, str(out)], capture_output=True, text=True)
        print(f"{name}: {hits} strings replaced | {r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr.strip()[-120:]}")
