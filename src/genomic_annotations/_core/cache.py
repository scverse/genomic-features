import pooch

PKG_CACHE_DIR = "genomic-annotations"


def retrieve_annotation(url):
    """Download and cache annotation file stored at url."""
    return pooch.retrieve(
        url=url,
        known_hash=None,
        path=pooch.os_cache(PKG_CACHE_DIR),
        progressbar=True,
    )
