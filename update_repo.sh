#!/usr/bin/env bash
# ==========================================================
# PhyloTracer 一键同步脚本（稳健版）
# - rebase 拉取，避免反复 merge commit
# - 任何一步失败直接退出并给出原因
# - 只有有改动才 commit
# ==========================================================

set -Eeuo pipefail

BRANCH="main"
REMOTE="origin"

log() { printf "🌿 %s\n" "$*"; }
warn() { printf "⚠️  %s\n" "$*" >&2; }
die() { printf "❌ %s\n" "$*" >&2; exit 1; }

# ---- precheck ----
[[ -d .git ]] || die "当前目录不是 Git 仓库，请进入项目根目录后再运行。"

# 确保在 main 分支（或你想要的分支）
current_branch="$(git branch --show-current)"
[[ "$current_branch" == "$BRANCH" ]] || die "当前分支是 $current_branch，请切换到 $BRANCH 再运行。"

# 如果上次 merge/rebase 中断，先退出
if [[ -d .git/rebase-apply || -d .git/rebase-merge ]]; then
  die "检测到正在进行 rebase。请先处理：git rebase --continue 或 git rebase --abort"
fi
if [[ -f .git/MERGE_HEAD ]]; then
  die "检测到正在进行 merge。请先处理：git status 查看冲突，然后 git merge --abort 或解决后提交"
fi

log "[PhyloTracer] Auto Git Sync Started..."

# ---- pull with rebase ----
log "🔄 Pulling latest changes from ${REMOTE}/${BRANCH} (rebase + autostash)..."
git fetch "$REMOTE" "$BRANCH"

# autostash：本地有未提交改动也能先暂存再 rebase
# 如果冲突会直接失败并退出（因为 set -e）
git pull --rebase --autostash "$REMOTE" "$BRANCH"

# ---- stage changes ----
log "📦 Adding all changed files..."
git add -A

# ---- commit if needed ----
if git diff --cached --quiet; then
  warn "No new changes to commit. (Index clean)"
else
  timestamp="$(date "+%Y-%m-%d %H:%M:%S")"
  commit_msg="AutoSync: update on $timestamp"
  log "📝 Committing changes: $commit_msg"
  git commit -m "$commit_msg"
fi

# ---- push (with retry if remote advanced) ----
log "🚀 Pushing to ${REMOTE}/${BRANCH}..."
if ! git push "$REMOTE" "$BRANCH"; then
  warn "Push failed. Remote may have advanced. Retrying once after rebase..."
  git fetch "$REMOTE" "$BRANCH"
  git pull --rebase --autostash "$REMOTE" "$BRANCH"
  git push "$REMOTE" "$BRANCH" || die "Push failed again. 请复制上面的 git 输出给我定位。"
fi

# ---- status ----
log "✅ Sync complete! Current status:"
git status -sb
log "--------------------------------------------------"
log "📁 Current branch: $(git branch --show-current)"
log "🌐 Remote: $(git remote get-url "$REMOTE")"
log "--------------------------------------------------"