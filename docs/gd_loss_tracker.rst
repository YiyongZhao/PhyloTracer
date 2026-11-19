GD_Loss_Tracker: counting modes and interpretation
=================================================

English
-------

Scope and state definitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The module maps each inferred gene duplication (GD) onto the species tree and
traces its two child lineages through descendant nodes. A numeric state of 2, 1
or 0 indicates how many of those child lineages are represented along the mapped
path. Each child contributes at most one, even if it contains multiple genes of
the same species. These states are not observed gene copy totals.

**Branch transitions** compare adjacent numeric path nodes:

* ``2_1`` (2→1): two represented GD child lineages decrease to one.
* ``1_0`` (1→0): the remaining represented lineage decreases to zero.
* ``2_0`` (2→0): two decrease directly to zero on one mapped path step, without
  an intervening state of 1. This does not prove simultaneous biological losses.
* ``2_2``, ``1_1`` and ``0_0`` describe unchanged states and are not loss
  transitions in the native node-loss summary.

Each decrease is one transition record. In particular, ``2_0`` contributes one
record rather than two lost-gene counts. For the illustrative path
``S20(2)->S19(1)->S18(1)->S17(0)``, the raw decreases are ``S19:2_1`` and
``S17:1_0``; S18 has no new decrease. ``S20(2)->S19(0)`` instead has one
``S19:2_0`` decrease. Node labels in these examples are illustrative and must be
checked against the numbered tree from the actual run.

**Terminal classification is a different object.** The column
``Species final loss type`` compares the two GD child subclades for one species:
``2_2`` means presence in both, ``2_1`` presence in one, and ``2_0`` absence from
both when that species is classifiable. Terminal ``2_0`` does not distinguish
2→0 from 2→1→0. There is no separate terminal ``1_0`` class. A species with
three genes in child A and one in child B is still terminal ``2_2``, not state 4.

Selecting a GD origin with --gd_node
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``--gd_node "SpeciesA,SpeciesB"`` selects GD events mapped **exactly to the
MRCA** of those two species and then traces the selected events through their
descendants. It does not select every GD below that ancestor, every path passing
through it, or every loss occurring there. Supply two names exactly as present
in the input species tree; quote the argument if a name contains spaces. Repeat
``--gd_node`` to select the union of several origin nodes. Generated S-node IDs
such as ``S20`` are not accepted by this option.

``--mrca_node`` is retained as a compatible legacy spelling for this module.
Hybrid_Tracer's separate ``--mrca_node`` option is unchanged.

Check the log's exact-MRCA message and the output ``GD burst node`` values.
Existing behavior skips invalid pairs; if no valid pair remains, no GD-origin
restriction is applied. This documentation change does not change that behavior.

This filter does not apply GD_Detector strict-event matching or any paper-specific
selection. A restricted rerun need not reproduce the event count of an independently
filtered dataset. To retain an already selected event set, filter its existing
CSV by ``GD burst node`` and keep its original Tree/GD identities.

Two counting modes
~~~~~~~~~~~~~~~~~~

``--node_count_mode parsimony`` (default)
    For each GD and transition type separately, collect species whose raw paths
    contain that type. Cover that affected set with maximal species-tree clades
    for which every descendant species is affected. Apply the internal-node
    support gates below and assign the type to each accepted clade root.
    The node summary counts each ``Tree ID + GD ID + Node + Type`` once.
    This is the implementation's clade-cover rule, not a general optimum over
    all possible evolutionary histories. Its assignments need not coincide with
    the first raw path decrease; ``GD loss pattern`` reflects those assignments.

``--node_count_mode accumulate``
    Keep each raw decrease on every descendant-species path and count every
    occurrence at its destination node. A shared ancestral decrease appearing
    on 19 exported species paths contributes 19 counts, not 19 independent
    biological loss events. The same Tree/GD can contribute repeatedly.

For one GD whose 19 descendant paths share a single accepted ancestral ``2_1``
assignment, parsimony contributes 1 at that assigned node; accumulate contributes
19 at the raw transition node if all 19 paths decrease there. The two nodes need
not always be the same because parsimony reassigns by clade cover.

Switching mode does not change GD detection, numeric ``GD loss path`` states,
A/B gene membership, or terminal classification with otherwise identical inputs.
It changes the effective loss assignments, ``GD loss pattern``, and node-summary
counts. ``GD loss node`` is derived separately from first raw decreases, not from
the parsimony assignments; do not assume it is identical to the nodes in the pattern.

Parsimony support gates
~~~~~~~~~~~~~~~~~~~~~~

* ``--parsimony_min_support_ratio`` (default 0.5) uses affected descendant species
  divided by **all descendant species of the candidate node**, not all species
  under the original GD node. It is not bootstrap or posterior branch support.
  The current clade-cover rule only selects fully affected clades, so this
  fraction is 1 for selected candidates; changing a cutoff within 0–1 does not
  make a partially affected clade eligible.
* ``--parsimony_min_support_species`` (default 2) sets the minimum affected-species
  count for internal candidate nodes. Terminal assignments are exempt.
* A candidate that fails either gate is skipped; it is not recursively re-expanded
  after rejection. Both parameters are ignored in accumulate mode.

Endpoint and missing-data policies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``--target_species SP`` keeps output rows and paths ending in SP; repeat for several
species. Transition inference is performed before this filter, so a restricted
export is not a fresh parsimony inference on only those terminal species.

By default, a species absent from the entire gene family receives terminal ``NA``
(internal reason ``missing_data``). With ``--include_unobserved_species``, it is
classified from GD-child presence instead. This flag changes terminal labels,
not GD detection or numeric path states. An ``NA`` row may still have a numeric
path ending in 0, and those numeric transitions can enter node counts in either
mode. ``NA`` is not automatically excluded from path-based calculations. Absence
in a family tree or both GD children is not independent proof of biological loss;
sampling, annotation and tree uncertainty remain relevant.

Outputs and cumulative histories
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``gd_loss.csv``
    One row per selected event/species record. Includes Tree/GD identities,
    ``GD burst node``, event-level ``GD loss node`` and ``GD loss pattern``, support,
    species, A/B gene lists and observed gene counts, terminal classification,
    path-expanded numeric columns, and ``GD loss path``. Event-level pattern
    strings deduplicate internal-node/type pairs and are repeated on species rows;
    summing their occurrences across rows does not recover either mode's counts.

``gd_loss_node_summary.tsv``
    Includes ``Node``, ``Mode``, ``GD event count`` (distinct Tree/GD origins at
    that node), ``Node loss count``, and either ``P2_0/P2_1/P1_0`` or
    ``C2_0/C2_1/C1_0``. Loss total is the sum of the three type columns. It does
    not include ``2_2`` retention or a complete inherited-state distribution.

``numed_sptree.nwk``
    The numbered species tree used to interpret the output labels. Generated
    numbering must not be assumed identical across different input trees.

A cumulative history of a fixed GD set is a **downstream summary**, not a third
``--node_count_mode`` or a native accumulate column. Select an explicit event set,
read each numeric path, and deduplicate by Tree/GD/target node after checking
repeated species paths agree. A state of 1 inherited from a parent remains 1,
without a new loss. If subdividing state 0 into direct 2→0 versus earlier
2→1→0 histories, inspect the preceding path states and document the policy for
``NA``. Such categories are histories, not adjacent-branch transitions. Their
total equals the selected GD count only when every selected GD has a consistent
state at the target node. The full gene-tree topology requires the original
rooted Newick tree; it cannot be reconstructed from the loss table alone.

Examples
~~~~~~~~

Use separate output directories for the two modes. The following example assumes
that SpeciesA and SpeciesB are valid names in the supplied tree:

.. code-block:: bash

    PhyloTracer GD_Loss_Tracker \
      --input_GF_list GF_ID2path.imap \
      --input_sps_tree sptree.nwk \
      --input_imap gene2sps.imap \
      --gd_node "SpeciesA,SpeciesB" \
      --node_count_mode parsimony \
      --output_dir loss_parsimony

    # Repeat the same inputs and GD-origin filter, changing only these options:
    # --node_count_mode accumulate --output_dir loss_accumulate

Keep each node-summary TSV next to its own CSV when using GD_Loss_Visualizer;
the visualizer can use that adjacent table to choose mode-specific counts.

中文译文
--------

状态和转换的含义
~~~~~~~~~~~~~~~~

每个 GD 有两个子谱系。路径中的 0/1/2 表示映射后有几个子谱系被表示，
不是实际基因条数；同一子谱系即使含多个同物种基因，也最多贡献 1。

* ``2_1``：相邻路径节点由 2 变为 1，即两个子谱系减少到一个。
* ``1_0``：相邻路径节点由 1 变为 0，即剩余的一个子谱系不再被表示。
* ``2_0``：相邻路径节点由 2 直接变为 0，中间没有状态 1。这是映射路径的
  一步下降，不证明两个生物学丢失同时发生。
* ``2_2``、``1_1``、``0_0`` 为状态不变，不计入原生节点丢失汇总。

每次下降计一条转换记录，``2_0`` 不按丢失两份基因加 2。
示例 ``S20(2)->S19(1)->S18(1)->S17(0)`` 在 S19 有 ``2_1``、
在 S17 有 ``1_0``，S18 没有新增下降。示例节点编号仅作说明。

终端列 ``Species final loss type`` 是另一种分类：``2_2`` 表示该物种在
两个 GD 子支都有基因，``2_1`` 表示只在一个子支有基因，``2_0`` 表示在
可分类的前提下两个子支均未观察到基因。终端 ``2_0`` 不区分直接 2→0
还是先 2→1 再 1→0，终端分类也没有独立的 ``1_0`` 类别。
例如 A 子支有 3 个、B 子支有 1 个基因，终端分类仍是 ``2_2``。

指定 GD 起源节点
~~~~~~~~~~~~~~~~

``--gd_node "物种A,物种B"`` 保留 GD 起源恰好映射在这两个物种 MRCA
上的事件，再追踪这些事件后代的路径。它不是筛选该节点的所有后代 GD、
经过该节点的所有路径，也不是筛选在该节点发生的 loss。
使用输入树中的准确物种名；名称含空格时对整个参数加引号。
可以重复指定以保留多个起源节点的并集，不接受直接输入 ``S20``。

本模块保留 ``--mrca_node`` 作为旧名兼容别名，Hybrid_Tracer 的同名参数不变。
检查日志中的 EXACT MRCA 信息和结果中的 ``GD burst node``。
现有实现会跳过无效物种对；若全部无效，就不限制 GD 起源。
本次仅改参数命名和说明，没有改变这一既有行为。

指定节点不会自动采用 GD_Detector strict 匹配或论文使用的额外筛选。
若要保持既有事件集，直接在其 CSV 中按 ``GD burst node`` 筛选，并保留原
Tree/GD 标识，比重跑后假定数量相同更可靠。

两种计数模式
~~~~~~~~~~~~

``parsimony``（默认，简约计数）
    对每个 GD、每种转换类型，找出原始路径出现过该转换的物种集合；
    用所有后代物种均属于该集合的最大完整类群覆盖它。
    通过内部节点支持门槛后，将转换分配到类群根；节点表按
    ``Tree ID + GD ID + Node + Type`` 去重，每个组合计一次。
    这是当前实现的类群覆盖规则，不代表穷举所有演化历史后的全局最优。
    分配节点不一定等于原始路径第一次下降的位置。

``accumulate``（路径累加计数）
    对每个后代物种的路径逐步统计原始下降。同一祖先下降若出现在
    19 条输出路径中，就计 19 次；这些不是 19 次独立的生物学丢失。

一个 GD 的 19 条路径若共享同一祖先下降，且满足简约分配条件，简约模式在
分配节点计 1；累加模式若 19 条路径均在同一节点下降，则在该原始节点计 19。
由于简约模式按类群重新定位，两模式的计数位置不保证总是一致。

切换模式不改变 GD 检测、数值路径、A/B 基因集合及同条件下的终端分类，
但改变有效转换分配、``GD loss pattern`` 和节点汇总。
``GD loss node`` 另由原始路径首次下降推得，不保证与简约 pattern 节点相同。

支持参数、缺失值及输出
~~~~~~~~~~~~~~~~~~~~~~

* ``--parsimony_min_support_ratio`` 默认 0.5，分母为候选节点的全部后代物种，
  不是整个 GD 起源类群，也不是基因树分支支持度。当前完整类群覆盖使入选
  候选的该比例为 1；在 0–1 内调整阈值不会使部分受影响类群入选。
* ``--parsimony_min_support_species`` 默认 2，只限制内部候选节点，终端豁免。
  不满足门槛的候选会被跳过，不会再拆回终端；两参数在 accumulate 中不使用。
* ``--target_species`` 在每个事件推断完成后筛选输出终点，不会在仅保留的
  物种集合上重新做一次简约推断。
* 默认将整个家族均未观察到的物种标为终端 ``NA``（内部原因 ``missing_data``）；
  ``--include_unobserved_species`` 改变终端分类策略，不改变数值路径。
  NA 行仍可能有终点为 0 的数值路径，相关下降仍可进入两模式的路径计数，
  不能假定 NA 已自动排除。基因树缺失和路径 0 本身不证明生物学丢失。

``gd_loss.csv`` 保存每个事件/物种的 A/B 基因、终端分类和路径；事件级
``GD loss pattern`` 是内部节点/类型的去重字符串，会重复写在物种行上，
不能跨行简单计数来代替节点表。
``gd_loss_node_summary.tsv`` 保存起源 GD 数、丢失转换总数及
``P2_0/P2_1/P1_0`` 或 ``C2_0/C2_1/C1_0``，不含 ``2_2`` 保留数。
``numed_sptree.nwk`` 保存对应节点编号，不同输入树的编号不保证可互换。

固定事件集在各节点的“累计历史”需要下游额外汇总，不是原生 accumulate
模式的直接输出：先固定事件集，再按 Tree/GD/目标节点去重，并检查重复物种
路径一致。继承的状态 1 不等于当前分支又发生一次丢失；把状态 0 细分为直接
2→0 或曾经过 1 的历史时，需读此前路径并声明 NA 政策。
只有每个事件在该目标节点都有一致状态时，分类和才等于固定事件集大小。
完整基因树拓扑仍需原始定根 Newick，不能仅靠 Tracker 表重建。
