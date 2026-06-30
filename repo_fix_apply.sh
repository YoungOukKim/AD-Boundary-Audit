#!/usr/bin/env bash
# =============================================================================
# AD-Boundary-Audit : 보충그림 번호를 원고(S1-S5)에 일치시키는 일괄 패치
# -----------------------------------------------------------------------------
# 원고 권위 매핑:
#   S1 pancellular | S2 GWAS | S3 global-correction-control(housekeeping)
#   S4 celltype-generalization | S5 pseudobulk | (EC: 무번호, Fig.3 보조)
# 기존 repo 문제: pseudobulk가 S3로 잘못 표기 + README에 S3 중복충돌
#                 (S3=pseudobulk vs S3 ctrl=housekeeping)
# -----------------------------------------------------------------------------
# 사용법: 저장소 루트(AD-Boundary-Audit/)에서
#         bash repo_fix_apply.sh
#         git diff --stat   # 검토
#         git commit -am "Align supplementary figure numbering with manuscript (S1-S5)"
#         git push
# (Windows는 Git Bash에서 실행)
# =============================================================================
set -euo pipefail
command -v git >/dev/null || { echo "git not found"; exit 1; }
[ -d .git ] || { echo "저장소 루트에서 실행하세요(.git 없음)"; exit 1; }

echo "[1/6] 스크립트 파일명 rename"
git mv R/figures/SuppFig_global_correction_control.R R/figures/SuppFigS3_global_correction_control.R
git mv R/figures/SuppFig_celltype_generalization.R   R/figures/SuppFigS4_celltype_generalization.R
git mv R/figures/SuppFigS3_pseudobulk.R              R/figures/SuppFigS5_pseudobulk.R

echo "[2/6] 스크립트 내부 출력경로/헤더/legend 갱신"
sed -i \
  -e 's/SuppFig_global_correction_control\.R/SuppFigS3_global_correction_control.R/g' \
  -e 's/SuppFig_global_correction_control\.png/SuppFigure_S3.png/g' \
  -e 's/SuppFig_global_correction_control\.tif/SuppFigure_S3.tif/g' \
  -e 's/SuppFig_global_correction_control\.{png,tif}/SuppFigure_S3.{png,tif}/g' \
  R/figures/SuppFigS3_global_correction_control.R
sed -i \
  -e 's/SuppFig_celltype_generalization\.R/SuppFigS4_celltype_generalization.R/g' \
  -e 's/SuppFig_celltype_generalization\.png/SuppFigure_S4.png/g' \
  -e 's/SuppFig_celltype_generalization\.tif/SuppFigure_S4.tif/g' \
  -e 's/SuppFig_celltype_generalization\.{png,tif}/SuppFigure_S4.{png,tif}/g' \
  R/figures/SuppFigS4_celltype_generalization.R
sed -i \
  -e 's/SuppFigS3_pseudobulk\.R/SuppFigS5_pseudobulk.R/g' \
  -e 's/SuppFigure_S3\.png/SuppFigure_S5.png/g' \
  -e 's/# Supplementary Fig\. S3\. Donor-level/# Supplementary Fig. S5. Donor-level/' \
  R/figures/SuppFigS5_pseudobulk.R

echo "[3/6] 커밋된 산출물 파일명 정리(이미지 내용 동일, 순서 주의)"
git mv output/figures/SuppFigure_S3.png                     output/figures/SuppFigure_S5.png
git mv output/figures/SuppFig_global_correction_control.png output/figures/SuppFigure_S3.png
git mv output/figures/SuppFig_global_correction_control.tif output/figures/SuppFigure_S3.tif
git mv output/figures/SuppFig_celltype_generalization.png   output/figures/SuppFigure_S4.png
git mv output/figures/SuppFig_celltype_generalization.tif   output/figures/SuppFigure_S4.tif
# 매뉴스크립트 reference 복사본(있으면) 번호 정정
if [ -f output/figures/_manuscript_reference/SuppFigure_S3_pseudobulk.png ]; then
  git mv output/figures/_manuscript_reference/SuppFigure_S3_pseudobulk.png \
         output/figures/_manuscript_reference/SuppFigure_S5_pseudobulk.png
fi

echo "[4/6] run_all.R 소스 경로 갱신"
sed -i \
  -e 's#R/figures/SuppFig_celltype_generalization\.R#R/figures/SuppFigS4_celltype_generalization.R#g' \
  -e 's#R/figures/SuppFig_global_correction_control\.R#R/figures/SuppFigS3_global_correction_control.R#g' \
  -e 's#R/figures/SuppFigS3_pseudobulk\.R#R/figures/SuppFigS5_pseudobulk.R#g' \
  run_all.R

echo "[5/6] README.md 보충그림 표 교체(S3 충돌 제거, S1-S5 정렬)"
python3 - << 'PY'
import re
p='README.md'; s=open(p,encoding='utf-8').read()
new=("| Supp. Fig. S1 | CPS-0.21 transition recovered across MTG cell types | `R/figures/SuppFigS1_pancellular.R` |\n"
"| Supp. Fig. S2 | AD-GWAS convergence on the CPS-0.207 boundary (microglia) | `R/figures/SuppFigS2_GWAS.R` |\n"
"| Supp. Fig. S3 | Housekeeping negative control for the global-expression correction | `R/figures/SuppFigS3_global_correction_control.R` |\n"
"| Supp. Fig. S4 | Audit generalizes across cell types (microglia panel) | `R/figures/SuppFigS4_celltype_generalization.R` |\n"
"| Supp. Fig. S5 | Donor-level (pseudobulk) testing is required | `R/figures/SuppFigS5_pseudobulk.R` |\n"
"| Supp. Fig. (EC) | Entorhinal non-replication (supports Fig. 3) | `R/figures/SuppFig_EC_nonreplication.R` |")
s2=re.sub(r"\| Supp\. Fig\. S1 \|.*?SuppFig_EC_nonreplication\.R` \|", new, s, flags=re.S)
assert s2!=s, "README supp-fig block 미발견 — 수동 확인 필요"
open(p,'w',encoding='utf-8').write(s2); print("  README.md OK")
PY

echo "[6/6] FIGURE_MAPPING.txt 교체(동봉 FIGURE_MAPPING.txt 사용 권장)"
echo "  -> 같이 받은 FIGURE_MAPPING.txt 로 덮어쓰세요 (cp /path/FIGURE_MAPPING.txt ./FIGURE_MAPPING.txt)"

echo ""
echo "완료. 'git diff --stat' 로 검토 후 commit/push 하세요."
