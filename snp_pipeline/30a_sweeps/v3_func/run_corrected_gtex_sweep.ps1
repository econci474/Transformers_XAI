# run_corrected_gtex_sweep.ps1
# ============================
# Drive the 63-run func-integration sweep on the corrected-GTEx
# AlphaGenome features (strict allow-list, 25-D, 2026-06-07 rerun).
#
# Grid (matches the 7 FM rows in the PRS master comparison table):
#   3 families × 7 func-integration modes × 3 seeds = 63 runs.
#
# Families:
#   - BMFM-SNP 2L MLP @ 101bp (chrom_hier, w=256, d=0.1, lr=0.001)
#   - NTv2     2L MLP @ 101bp (chrom_hier, w=512, d=0.3, lr=0.0003)
#   - BMFM-SNP XGB    @ 1001bp (global_attn, xgb_lr=0.05, depth=4, n_est=20)
#
# All runs to:
#   wandb project = snp_gtex_corrected_func_sweep
#   wandb group   = strict_v1_layerhippo
#   --ag-parent   = outputs/corrected_gtex/alphagenome_features
#   --output-root = outputs/corrected_gtex/diff_attn_v3
#
# Usage:
#   pwsh Transformers_XAI/snp_pipeline/30a_sweeps/v3_func/run_corrected_gtex_sweep.ps1
#   pwsh ... -SkipNone    # skip the 9 `none` control runs (assumes original values)
#   pwsh ... -SeedSubset 0,1   # only seeds 0 + 1
#   pwsh ... -Families bmfm_mlp,ntv2_mlp   # subset of the 3 families

param(
    [switch] $SkipNone,
    [int[]]  $SeedSubset = @(0, 1, 2),
    [string[]] $Families = @("bmfm_mlp", "ntv2_mlp", "bmfm_xgb")
)

$REPO = "C:/Users/elena/iCloudDrive/Desktop/ACS_MPhil/Thesis/git"
$PY   = "C:/Users/elena/miniforge3/envs/snp/python.exe"
$AG   = "$REPO/Transformers_XAI/snp_pipeline/outputs/corrected_gtex/alphagenome_features"
$OUT  = "$REPO/Transformers_XAI/snp_pipeline/outputs/corrected_gtex/diff_attn_v3"
$WB_PROJECT = "snp_gtex_corrected_func_sweep"
$WB_GROUP   = "strict_v1_layerhippo"

$MODES = @(
    "none",
    "per_snp_zscored_25",
    "per_snp_summaries_12",
    "attn_bias_single",
    "attn_bias_per_modality",
    "value_mult",
    "attn_bias_per_modality_plus_value_mult",
    "attn_bias_per_modality_plus_per_snp_summaries_12",
    "attn_bias_single_plus_func_imp_abs"
)
if ($SkipNone) { $MODES = $MODES | Where-Object { $_ -ne "none" } }

$FAMILY_CFGS = @{
    "bmfm_mlp" = @{
        args = @("--seq-length", "101bp", "--model", "bmfm_snp",
                  "--head", "mlp2", "--aggregation", "chrom_hier",
                  "--mlp-width", "256", "--mlp-dropout", "0.1",
                  "--lr", "0.001")
        label = "BMFM-SNP 2L MLP @ 101bp"
    }
    "ntv2_mlp" = @{
        args = @("--seq-length", "101bp", "--model", "ntv2",
                  "--head", "mlp2", "--aggregation", "chrom_hier",
                  "--mlp-width", "512", "--mlp-dropout", "0.3",
                  "--lr", "0.0003")
        label = "NTv2 2L MLP @ 101bp"
    }
    "bmfm_xgb" = @{
        args = @("--seq-length", "1001bp", "--model", "bmfm_snp",
                  "--head", "xgb", "--aggregation", "global_attn",
                  "--xgb-lr", "0.05", "--xgb-max-depth", "4",
                  "--xgb-n-estimators", "20")
        label = "BMFM-SNP XGB @ 1001bp"
    }
}

$total = $Families.Count * $MODES.Count * $SeedSubset.Count
$i = 0
$t0 = Get-Date
Write-Host "=== Corrected-GTEx func-integration sweep ===" -ForegroundColor Cyan
Write-Host "Total planned: $total runs ($($Families.Count) families x $($MODES.Count) modes x $($SeedSubset.Count) seeds)"
Write-Host "WandB project: $WB_PROJECT  /  group: $WB_GROUP"
Write-Host "Output root:   $OUT"
Write-Host "AG parent:     $AG"
Write-Host ""

function _fmtElapsed($ts) {
    return ('{0:00}h{1:00}m{2:00}s' -f $ts.Hours, $ts.Minutes, $ts.Seconds)
}

foreach ($seed in $SeedSubset) {
    foreach ($mode in $MODES) {
        foreach ($fam in $Families) {
            $i++
            $cfg = $FAMILY_CFGS[$fam]
            $elapsed = _fmtElapsed((Get-Date) - $t0)
            Write-Host "[$i/$total] $($cfg.label)  mode=$mode  seed=$seed  (elapsed: $elapsed)" -ForegroundColor Yellow
            $cmd = @($PY, "-X", "utf8",
                      "$REPO/Transformers_XAI/snp_pipeline/30v3_train_diff_attention_func.py",
                      "--ag-parent", $AG,
                      "--output-root", $OUT,
                      "--func-integration-mode", $mode,
                      "--seed", "$seed",
                      "--wandb-project", $WB_PROJECT,
                      "--wandb-group", $WB_GROUP) + $cfg.args
            & $cmd[0] $cmd[1..($cmd.Length - 1)]
            if ($LASTEXITCODE -ne 0) {
                Write-Host "  [FAIL] exit code $LASTEXITCODE -- continuing" -ForegroundColor Red
            }
        }
    }
}

$total_time = _fmtElapsed((Get-Date) - $t0)
Write-Host ""
Write-Host "=== Sweep complete: $i runs in $total_time ===" -ForegroundColor Green
