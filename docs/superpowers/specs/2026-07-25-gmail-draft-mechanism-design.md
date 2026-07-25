# Gmail draft mechanism — design

**Date:** 2026-07-25
**Status:** ON HOLD — build-vs-adopt not yet settled
**Decision:** Approach A — Gmail API with `gmail.compose`, wrapped in a tool that has no send path.

> **Hold, added after writing.** Ben raised that a prebuilt MCP server or CLI
> probably already does this, and he is very likely right. Purpose-built Gmail
> MCP servers exist, several with OAuth handled for you and a draft tool
> exposed directly. **If a maintained one checks out, this spec should be
> abandoned rather than implemented** — writing a bespoke OAuth client to
> duplicate a maintained package is the wrong trade, and every line of it is
> credential-handling code that would need reviewing.
>
> The research was blocked by a classifier outage. Do it before writing code:
>
> 1. `mcp__mcp-registry__search_mcp_registry` for gmail / email / workspace.
> 2. Search for maintained Gmail MCP servers; check last commit, whether
>    `drafts.create` is exposed, and **which scopes they request** — many ask
>    for `gmail.modify` or full access, which is far broader than needed here
>    and would defeat the point.
> 3. CLI fallbacks worth pricing: `himalaya` (Rust, IMAP/SMTP), `gam` /
>    `GAMADV-XTD3` (admin-oriented, heavy).
>
> **What does not change either way:** most self-hosted options still need the
> same Google Cloud project and OAuth desktop client — roughly the same 20
> minutes of console clicks in the Setup section below. Hosted options (Zapier,
> Composio and similar) avoid that but put a third party in the path of your
> mail, which is a materially different trust decision and should be taken
> deliberately, not by default.
>
> The two findings below survive regardless of what gets adopted: the scope
> question, and the verification requirement. Apply both to any adopted tool —
> in particular, **check the scopes it requests rather than assuming they are
> minimal.**

## Problem

Creating a Gmail draft from an agent session currently has no working mechanism.
Three attempts failed today, each in a different way:

- **No Gmail connector.** Tool search returns Slack, Sanity and browser tools only.
- **Claude in Chrome not connected.** `list_connected_browsers` returns empty.
- **`open -a "Google Chrome" <compose-url>`** reported success and did nothing.
  `open` only dispatches; it confirms neither launch nor navigation. Chrome was
  not running afterwards.

That last one is the instructive failure. The command's exit code was reported
as success when nothing had been verified — the same mistake as checking that a
link *exists* rather than where it sits on the page. **Verification of the
outcome, not the call, is a first-class requirement of this design.**

Three letters are queued behind this (White, SimpleTES authors, Kim & Pilanci),
so a reusable mechanism is worth the setup.

## The scope finding, and a hard precondition

The requested constraint was "drafts only, no send scope." **I do not believe
that constraint is expressible in the Gmail API.** `gmail.compose` is, to my
knowledge, documented as managing drafts *and sending*.

**PRECONDITION — verify before writing any code.** Fetch
`https://developers.google.com/gmail/api/auth/scopes` and confirm the exact
capabilities of `gmail.compose`. Then:

- **If `gmail.compose` permits sending** (expected): proceed as designed, and
  state the security property honestly as "the tool has no send path," not
  "Google prevents sending."
- **If a narrower draft-only scope exists**: use it instead, and this document's
  security section is simply stronger than written.

This was blocked today by a classifier outage. It is not optional.

### Why `gmail.compose` is still the narrowest sensible credential

It **cannot read the mailbox.** A token that writes drafts but cannot list,
search, or read a single existing message is genuinely narrow, and revocable in
one click from the Google account page.

The rejected alternative makes this concrete: IMAP `APPEND` *structurally*
cannot send, which sounds better — but it requires an app password, and an app
password grants full IMAP **and** SMTP. Narrow operation, broadest possible
credential. Rejected on least-privilege.

## Architecture

One module, one job.

```
~/.local/share/claude-gmail-draft/
    gmail_draft.py         the tool (executable, own venv shebang)
~/.local/bin/
    gmail-draft            symlink -> ../share/claude-gmail-draft/gmail_draft.py
~/.config/claude-gmail-draft/
    client_secret.json     OAuth desktop client        (chmod 600)
    token.json             refresh token, created on first run (chmod 600)
```

**Deliberately outside any git repository.** The token cannot be committed by
accident because it does not live where commits happen. This is not a
`.gitignore` promise; it is a location.

Dependencies (`google-auth-oauthlib`, `google-api-python-client`) install into a
dedicated venv at `~/.local/share/claude-gmail-draft/.venv`, not into the Erdős
project venv — this tool has nothing to do with that project and should not be
able to break it.

### Interface

`gmail-draft` is a symlink to `gmail_draft.py`; the two names refer to one file.

```
gmail-draft --to ADDR --subject TEXT --body-file PATH [--thread-id ID] [--dry-run]
```

- `--dry-run` builds the MIME and prints it. No network, no auth. This is the
  default posture for any first use of a new letter.
- `--thread-id` makes the draft a reply in an existing thread, which the White
  letter needs.
- Exit non-zero on any failure. Never print success without having verified.

### Flow

1. Load token; if absent or expired, run `InstalledAppFlow` against a loopback
   listener. Browser opens once, consent granted, token written.
2. Build a MIME message. **UTF-8 throughout** — the letters contain `µ`, `≥`,
   and em-dashes, and mangled encoding is a realistic failure here.
3. `drafts.create`.
4. **Verify:** `drafts.get` on the returned id, and assert the recipient and
   subject round-tripped. Print the draft id and its Gmail URL.

Step 4 is the point of the design. Today's failures were all "the call
returned, so I said it worked."

## Security properties

Stated precisely, because the honest version is weaker than the one requested:

| Property | Guaranteed by |
| --- | --- |
| Cannot read your mail | Google — `gmail.compose` has no read capability |
| Cannot send | **The tool only.** No send code path exists in it |
| Token not committed | Location outside any repo, plus `chmod 600` |
| Revocable | One click, Google account permissions page |

The send guarantee is behavioural, not enforced. If that is unacceptable, the
Claude in Chrome extension is the only option with no standing credential — at
the cost of needing Chrome connected and offering no unattended operation.

## Error handling

- **No token, non-interactive session:** fail with the exact command to run
  interactively. Never hang waiting on a browser that cannot open.
- **Token expired/revoked:** detect, delete the stale token, tell the user to
  re-run once interactively.
- **API error:** surface Google's message verbatim. Do not retry blindly — a
  403 means a scope or consent problem that retrying cannot fix.
- **Verification mismatch:** treat as failure, exit non-zero, print both what
  was sent and what came back.

## Testing

- `--dry-run` on the existing White letter; diff the MIME body against
  `communications/email_to_ethan_white_v6.txt`.
- **Unicode round-trip:** draft containing `µ ≥ 0.3803954` and an em-dash,
  verified via `drafts.get`, asserting exact equality. This is the most likely
  real-world break.
- **Self-test:** create a draft addressed to Ben, verify, then delete it.
  Leaves no residue.
- **No send path:** `grep -c "send" gmail_draft.py` as a crude but real check
  that the capability was never added later.

## Setup — who does what

**Ben (~20 min, cannot be delegated):** create a Google Cloud project, enable
the Gmail API, create an OAuth **Desktop** client, download the client secret to
`~/.config/claude-gmail-draft/client_secret.json`, add himself as a test user on
the consent screen. Then run the tool once interactively to grant consent.

**Me:** everything else — the tool, tests, and the first real draft.

## Out of scope

Sending. Reading mail. Attachments. HTML bodies (plain text only; these are
letters to mathematicians). Multiple accounts.
