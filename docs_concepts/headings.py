from pathlib import Path
from headings_util import print_md_toc, split_markdown_by_heading_ranges, print_md_section_line_counts, first_md_starting


if __name__ == "__main__":
    script_dir = Path(__file__).parent

    # You can now enter the start of the filename here:
    filename_start = "briefing"
    md_path = first_md_starting(script_dir, filename_start)

    print = 0

    if print:
        print_md_section_line_counts(md_path, level=2)
    else:
        split_markdown_by_heading_ranges(md_path, [6,8,10,11,12,13,14,15,16,17], heading_level=2)


