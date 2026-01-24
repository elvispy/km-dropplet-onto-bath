# .context Protocol (Authoritative)
- .context is agent-only working memory. Do not quote or dump it to the user unless explicitly asked.
- Every task must have exactly one task file: .context/tasks/<task_id>.md
- Work must be traceable: all non-trivial actions update the task file (plan/log/handoff).
- Inter-agent communication happens via .context/INBOX.md and task-file “Messages” sections.
- Append-only rule: never delete other agents’ notes; prefer appending with timestamps.
- Locking rule: one active owner per task at a time (lease). If you need to take over, record the claim in TASK_INDEX + task file.
