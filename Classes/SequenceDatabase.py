class SequenceDatabase(dict):
    def __getitem__(self, key):
        if key not in self:
            raise KeyError(f"ERROR! Missing sequence: {key}")
        return super().__getitem__(key)
