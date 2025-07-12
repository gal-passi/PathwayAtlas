from definitions import DIRS_TO_CREATE

def setup_directories():
    for path in DIRS_TO_CREATE:
        os.makedirs(path, exist_ok=True)

if __name__ == "__main__":
    setup_directories()