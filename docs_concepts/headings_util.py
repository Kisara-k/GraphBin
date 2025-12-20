def print_md_toc(md_file_path: str, level: int = 3):
    try:
        with open(md_file_path, "r", encoding="utf-8") as f:
            for line_number, line in enumerate(f, start=1):
                stripped = line.lstrip()
                if stripped.startswith("#"):
                    heading_level = len(stripped) - len(stripped.lstrip("#"))
                    if heading_level <= level:
                        heading_text = stripped[heading_level:].strip()
                        indent = "  " * (heading_level - 1)
                        print(f"{line_number:03d} {indent}- {heading_text}")
    except FileNotFoundError:
        print(f"Markdown file not found: {md_file_path}")

def print_md_section_line_counts(md_file_path: str, level: int):
    """
    Prints the number of lines under each heading at EXACTLY the given level.
    Sections only end when encountering another heading at the SAME level.
    Higher-level (parent) headings don't end sections, they're just part of the text.
    Lower-level (child) headings also don't end sections, they're part of the section content.
    """

    sections = []
    current_section = None

    with open(md_file_path, "r", encoding="utf-8") as f:
        for line in f:
            stripped = line.lstrip()

            # Check if this line is a heading
            if stripped.startswith("#"):
                heading_level = len(stripped) - len(stripped.lstrip("#"))

                # If we encounter another heading at the SAME level,
                # end the current section and start a new one
                if heading_level == level:
                    if current_section is not None:
                        sections.append(current_section)
                    
                    # Start new section
                    heading_text = stripped[heading_level:].strip()
                    current_section = {
                        "title": heading_text,
                        "lines": 0
                    }
                # Don't do anything for other heading levels - 
                # they're just part of the text content for our counting

                continue

            # Count lines only if we are inside a target-level section
            if current_section is not None:
                current_section["lines"] += 1

        # Close last section if file ends
        if current_section:
            sections.append(current_section)

    # Output
    idx = 1
    for s in sections:
        print(f"{s['lines']:4d} | {idx:02d} | {s['title']}")
        idx += 1

from pathlib import Path

from pathlib import Path
import re

def split_markdown_by_heading_ranges(
    md_file_path: str,
    split_points: list[int],
    heading_level: int = 1,
):
    """
    Split a markdown file into multiple files based on heading index ranges.
    Each output file is named:
      NN-heading-title.md
    where NN is the split point (zero-padded).
    """

    if type(split_points) is tuple and len(split_points) == 2:
        a = split_points[0]
        b = split_points[1]
        split_points = list(range(a, b + 2))

    md_path = Path(md_file_path)
    output_dir = md_path.with_suffix("")
    output_dir.mkdir(exist_ok=True)

    heading_prefix = "#" * heading_level + " "

    lines = md_path.read_text(encoding="utf-8").splitlines(keepends=True)

    # Collect heading line indices and titles
    headings = []
    for i, line in enumerate(lines):
        if line.startswith(heading_prefix):
            title = line[len(heading_prefix):].strip()
            headings.append((i, title))

    if not headings:
        return

    total_headings = len(headings)

    # Build ranges: (start_heading_idx, end_heading_idx, label_number)
    ranges = []
    prev = 1

    for point in split_points:
        ranges.append((prev, point - 1, prev))
        prev = point

    ranges.append((prev, total_headings, prev))

    for start_h, end_h, label in ranges:
        if start_h > total_headings:
            break

        start_h = max(1, start_h)
        end_h = min(end_h, total_headings)

        start_line, title = headings[start_h - 1]
        end_line = (
            headings[end_h][0]
            if end_h < total_headings
            else len(lines)
        )

        content = lines[start_line:end_line]

        # Sanitize heading title for filename
        safe_title = re.sub(r"[^\w\-]+", "_", title).strip("_")

        filename = f"{label:02d}-{safe_title}.md"
        (output_dir / filename).write_text(
            "".join(content),
            encoding="utf-8"
        )



def first_md_starting(directory: Path, start: str) -> Path | None:
    """Return the first .md file in the directory that starts with `start` (case-insensitive)."""
    start_lower = start.lower()
    for file in directory.glob("*.md"):
        if file.name.lower().startswith(start_lower):
            return file
    return None