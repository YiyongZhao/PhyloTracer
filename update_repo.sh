#!/bin/bash
# ==========================================================
# PhyloTracer 一键同步脚本
# 自动执行 git add → commit → push → 状态报告
# ==========================================================

echo "🌿 [PhyloTracer] Auto Git Sync Started..."

# 检查是否在 Git 仓库中
if [ ! -d .git ]; then
  echo "❌ Error: 当前目录不是一个 Git 仓库。请进入项目根目录后再运行。"
  exit 1
fi

# 拉取远程更新，防止冲突
echo "🔄 Pulling latest changes from origin/main..."
git fetch origin main
git merge origin/main --no-edit

# 添加所有变动（新增/删除/修改）
echo "📦 Adding all changed files..."
git add -A

# 自动生成提交信息
timestamp=$(date "+%Y-%m-%d %H:%M:%S")
commit_msg="AutoSync: update on $timestamp"

# 提交变动
echo "📝 Committing changes..."
git commit -m "$commit_msg" 2>/dev/null || echo "⚠️ No new changes to commit."

# 推送到 GitHub
echo "🚀 Pushing to origin/main..."
git push origin main

# 显示状态
echo "✅ Sync complete! Current status:"
git status -s
echo "--------------------------------------------------"
echo "📁 Current branch: $(git branch --show-current)"
echo "🌐 Remote: $(git remote get-url origin)"
echo "--------------------------------------------------"

