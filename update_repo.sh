#!/usr/bin/env bash
set -Eeuo pipefail

BRANCH="main"
REMOTE="origin"

echo "🌿 [PhyloTracer] Auto Git Sync Started..."

if [ ! -d .git ]; then
  echo "❌ Error: 当前目录不是 Git 仓库。"
  exit 1
fi

echo "🔄 Pulling latest changes..."
git fetch $REMOTE $BRANCH
git pull --rebase --autostash $REMOTE $BRANCH

echo "📦 Adding files..."
git add -A

if git diff --cached --quiet; then
  echo "⚠️ No changes to commit."
else
  timestamp=$(date "+%Y-%m-%d %H:%M:%S")
  commit_msg="AutoSync: update on $timestamp"
  echo "📝 Committing..."
  git commit -m "$commit_msg"
fi

echo "🚀 Pushing to $REMOTE/$BRANCH..."

if git push $REMOTE $BRANCH; then
  echo "🟢 PUSH SUCCESS"
  
  # 显示远程 commit
  remote_commit=$(git ls-remote $REMOTE $BRANCH | cut -f1)
  echo "🌐 Remote HEAD: $remote_commit"
  
else
  echo "🔴 PUSH FAILED"
  echo "⚠️ 可能原因："
  echo "   - 远程有更新未 pull"
  echo "   - 权限问题"
  echo "   - 网络问题"
  exit 1
fi

echo "--------------------------------------------------"
git status -sb
echo "--------------------------------------------------"