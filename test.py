import subprocess

def test_tool(path, args):
    try:
        subprocess.run([path] + args, check=True)
        print(f"{path} działa poprawnie.")
    except Exception as e:
        print(f"Błąd: {path} nie działa. Szczegóły: {e}")

# Ścieżki do narzędzi
P2RANK_PATH = "/usr/local/bin/prank"
VINA_PATH = "/usr/local/bin/vina_1.2.5_linux_x86_64"
OBABEL_PATH = "/usr/bin/obabel"

# Testowanie narzędzi
test_tool(P2RANK_PATH, ['--help'])
test_tool(VINA_PATH, ['--version'])
test_tool(OBABEL_PATH, ['--version'])
