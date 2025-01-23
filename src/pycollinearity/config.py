import argparse


def add_options(parser):
    parser.add_argument('--ref', type=str, required=True, help='reference file to index')
    parser.add_argument('--background', action='store_true', default=False,
                        help='use this flag if the index is of background species (whose reads will be discarded)')
    parser.add_argument('--minknow-addr', type=str, default='localhost:8000', help='minknow address')
    parser.add_argument('--basecaller', type=str, default='dorado', help='basecaller', choices=['dorago', 'guppy', ])
    parser.add_argument('--basecaller-addr', type=str, default='ipc:///tmp/.guppy/5555', help='basecaller address')
    parser.add_argument('--device', type=str, required=True, help='Nanopore device name')


class Config:
    def __init__(self):
        parser = argparse.ArgumentParser()
        add_options(parser)
        self.args = parser.parse_args()

    def __getattr__(self, attr):
        return getattr(self.args, attr)


if __name__ == '__main__':
    config = Config()
    print(config.ref)
    print(config.basecaller)
    print(config.basecaller_addr)
