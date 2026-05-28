#!/usr/bin/env python3
"""
Upload Three-phase methodology markdown docs to Feishu cloud documents.
Converts markdown → Feishu Docx blocks via Open API.
"""
import json, re, sys, io
import urllib.request
import urllib.error

# Fix Windows GBK encoding issues
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# ── Config ──────────────────────────────────────────────
APP_ID = "cli_aa838bc6a5785cc4"
APP_SECRET = "gRhU1cBBf9j6Kg10gVGTseOrrnwoTZVw"
TOKEN_URL = "https://open.feishu.cn/open-apis/auth/v3/tenant_access_token/internal"
DOC_CREATE_URL = "https://open.feishu.cn/open-apis/docx/v1/documents"
BLOCK_CHILDREN_URL = "https://open.feishu.cn/open-apis/docx/v1/documents/{document_id}/blocks/{block_id}/children"
FOLDER_TOKEN = ""  # set to a folder token to create docs inside a folder

DOCS = [
    {
        "file": "01_Flory-Huggins三元相图计算方法论.md",
        "title": "Flory-Huggins 三元相图计算方法论",
        "desc": "OPV三元体系相分离计算——FH自由能、Spinodal、Critical Point、Binodal全流程"
    },
    {
        "file": "02_HSP_相图集成溶剂选择框架.md",
        "title": "HSP + 三元相图集成溶剂选择框架",
        "desc": "三层筛选逻辑——从HSP热力学匹配到相图动力学再到蒸发补偿的完整决策框架"
    },
    {
        "file": "03_FH参数优化协议.md",
        "title": "Flory-Huggins 参数优化协议",
        "desc": "三维网格扫描策略、评分函数设计、参数调优规律与新体系参数化流程"
    },
]

# ── Auth ────────────────────────────────────────────────
def get_token():
    data = json.dumps({"app_id": APP_ID, "app_secret": APP_SECRET}).encode()
    req = urllib.request.Request(TOKEN_URL, data=data, headers={"Content-Type": "application/json"})
    resp = json.loads(urllib.request.urlopen(req).read())
    if resp.get("code") != 0:
        raise RuntimeError(f"Token error: {resp}")
    return resp["tenant_access_token"]

TOKEN = get_token()
print("Token obtained OK")

# ── API helpers ─────────────────────────────────────────
def api_post(url, body, method="POST"):
    data = json.dumps(body).encode() if body else None
    req = urllib.request.Request(url, data=data, method=method)
    req.add_header("Authorization", f"Bearer {TOKEN}")
    req.add_header("Content-Type", "application/json")
    try:
        resp = json.loads(urllib.request.urlopen(req).read())
        if resp.get("code") != 0:
            print(f"  API ERROR: {resp.get('code')} - {resp.get('msg')}")
        return resp
    except urllib.error.HTTPError as e:
        print(f"  HTTP ERROR {e.code}: {e.read().decode()[:300]}")
        return None

def create_doc(title, folder_token=""):
    """Create a new Feishu Docx document. Returns document_id."""
    body = {"title": title}
    if folder_token:
        body["folder_token"] = folder_token
    resp = api_post(DOC_CREATE_URL, body)
    if resp and resp.get("code") == 0:
        doc_id = resp["data"]["document"]["document_id"]
        return doc_id
    return None

def add_blocks(document_id, parent_block_id, blocks):
    """Add blocks as children of parent_block_id."""
    body = {"children": blocks}
    url = BLOCK_CHILDREN_URL.format(document_id=document_id, block_id=parent_block_id)
    return api_post(url, body)

# ── Markdown → Feishu Blocks ────────────────────────────
def text_element(text, bold=False, italic=False, code=False):
    """Create a Feishu text element with optional styles."""
    el = {"text_run": {"content": text}}
    style = {}
    if bold: style["bold"] = True
    if italic: style["italic"] = True
    if code: style["inline_code"] = True
    if style:
        el["text_run"]["text_element_style"] = style
    return el

def parse_inline(line):
    """Parse inline markdown (bold, code, italic) → list of text_elements."""
    elements = []
    # Handle **bold**
    parts = re.split(r'(\*\*.*?\*\*)', line)
    for part in parts:
        if part.startswith('**') and part.endswith('**'):
            elements.append(text_element(part[2:-2], bold=True))
        elif part.startswith('*') and part.endswith('*') and not part.startswith('**'):
            elements.append(text_element(part[1:-1], italic=True))
        else:
            # Handle inline `code`
            subparts = re.split(r'(`[^`]+`)', part)
            for sp in subparts:
                if sp.startswith('`') and sp.endswith('`'):
                    elements.append(text_element(sp[1:-1], code=True))
                else:
                    elements.append(text_element(sp))
    return [e for e in elements if e["text_run"]["content"]]

def remove_heading_marks(text):
    """Remove ## markers and bold marks from heading text."""
    t = re.sub(r'^#+\s*', '', text)
    t = re.sub(r'\*\*(.*?)\*\*', r'\1', t)
    return t

def make_block(block_type, elements=None):
    """Helper to create a block dict."""
    b = {"block_type": block_type}
    if block_type in (2, 3, 4, 5, 6):
        key = {2:"text", 3:"heading1", 4:"heading2", 5:"heading3", 6:"heading4"}[block_type]
        b[key] = {"elements": elements or [], "style": {}}
    elif block_type == 14:  # code block
        b["code"] = {"elements": elements or [], "style": {"language": 1}}
    elif block_type == 22:  # divider
        b["divider"] = {}
    return b

def md_to_blocks(md_text):
    """Convert markdown text to a list of Feishu block dicts."""
    blocks = []
    lines = md_text.split('\n')
    i = 0

    while i < len(lines):
        line = lines[i]

        # Skip empty lines
        if not line.strip():
            i += 1
            continue

        # Code blocks: ``` ... ```
        if line.strip().startswith('```'):
            i += 1
            code_lines = []
            while i < len(lines) and not lines[i].strip().startswith('```'):
                code_lines.append(lines[i])
                i += 1
            i += 1  # skip closing ```
            code_text = '\n'.join(code_lines)
            if code_text.strip():
                blocks.append(make_block(14, [text_element(code_text)]))
            continue

        # --- divider
        if re.match(r'^---+$', line.strip()):
            blocks.append(make_block(22))
            i += 1
            continue

        # ## Headings
        if line.startswith('#### '):
            blocks.append(make_block(6, parse_inline(remove_heading_marks(line))))
        elif line.startswith('### '):
            blocks.append(make_block(5, parse_inline(remove_heading_marks(line))))
        elif line.startswith('## '):
            blocks.append(make_block(4, parse_inline(remove_heading_marks(line))))
        elif line.startswith('# '):
            blocks.append(make_block(3, parse_inline(remove_heading_marks(line))))
        # Table rows: leave as formatted text
        elif line.strip().startswith('|'):
            blocks.append(make_block(2, [text_element(line.strip())]))
        # Bullet list → text block with bullet prefix
        elif re.match(r'^[\-\*]\s+', line.strip()):
            content = re.sub(r'^[\-\*]\s+', '• ', line.strip())
            blocks.append(make_block(2, parse_inline(content)))
        # Ordered list → text block with number prefix preserved
        elif re.match(r'^\d+\.\s+', line.strip()):
            blocks.append(make_block(2, parse_inline(line.strip())))
        # Regular paragraph
        else:
            clean = line.strip()
            # Skip markdown reference links like [text][id]
            if clean.startswith('[') and ']:' in clean:
                i += 1
                continue
            blocks.append(make_block(2, parse_inline(clean)))

        i += 1

    return blocks

# ── Main ────────────────────────────────────────────────
def main():
    import os
    script_dir = os.path.dirname(os.path.abspath(__file__))
    created_docs = []

    for doc_info in DOCS:
        title = doc_info["title"]
        filepath = os.path.join(script_dir, doc_info["file"])

        if not os.path.exists(filepath):
            print(f"SKIP {title} — file not found: {filepath}")
            continue

        with open(filepath, 'r', encoding='utf-8') as f:
            md_text = f.read()

        print(f"\nProcessing: {title}")
        print(f"  Markdown: {len(md_text)} chars")

        # Create document
        doc_id = create_doc(title, FOLDER_TOKEN)
        if not doc_id:
            print(f"  FAILED to create document")
            continue
        print(f"  Doc ID: {doc_id}")

        # Get the document's root block
        # The document_id itself is the parent for top-level blocks
        # First, get the document info to find the page block
        doc_info_url = f"https://open.feishu.cn/open-apis/docx/v1/documents/{doc_id}"
        req = urllib.request.Request(doc_info_url)
        req.add_header("Authorization", f"Bearer {TOKEN}")
        req.add_header("Content-Type", "application/json")
        doc_resp = json.loads(urllib.request.urlopen(req).read())
        if doc_resp.get("code") != 0:
            print(f"  FAILED to get doc info: {doc_resp}")
            continue

        # The document has a root block_id (usually same as document_id)
        page_block_id = doc_resp["data"]["document"].get("block_id", doc_id)
        print(f"  Page block: {page_block_id}")

        # Convert markdown to blocks
        blocks = md_to_blocks(md_text)
        print(f"  Blocks: {len(blocks)}")

        # Add blocks in batches of 50 (API limit)
        batch_size = 50
        for start in range(0, len(blocks), batch_size):
            batch = blocks[start:start + batch_size]
            resp = add_blocks(doc_id, page_block_id, batch)
            if resp and resp.get("code") == 0:
                print(f"  Batch {start//batch_size + 1}/{(len(blocks)-1)//batch_size + 1} OK ({len(batch)} blocks)")
            else:
                print(f"  Batch {start//batch_size + 1} FAILED")

        # URL
        url = f"https://{doc_resp['data']['document']['url'] if 'url' in doc_resp.get('data', {}).get('document', {}) else 'bytedance.feishu.cn/docx/' + doc_id}"
        created_docs.append({"title": title, "doc_id": doc_id, "url": f"https://bytedance.feishu.cn/docx/{doc_id}"})
        print(f"  URL: https://bytedance.feishu.cn/docx/{doc_id}")

    # ── Summary ──
    print("\n" + "=" * 60)
    print("UPLOAD SUMMARY")
    print("=" * 60)
    for d in created_docs:
        print(f"  {d['title']}")
        print(f"  → {d['url']}")
        print()

if __name__ == "__main__":
    main()
