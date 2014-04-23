
_command_classes = { }


class Meta(type):
    def __new__(cls, name, bases, attrs):
        newcls = super(Meta, cls).__new__(cls, name, bases, attrs)
        if not attrs.get('_skip', False):
            _command_classes[name] = newcls
        return newcls


class Command(object):
    __metaclass__ = Meta
    _skip = True

    def configure_parser(self, parser):
        pass


def create_parser(description=None):
    import argparse
    parser = argparse.ArgumentParser(description=description)
    subcmd = parser.add_subparsers()

    for name, cmdcls in _command_classes.iteritems():
        cmd = cmdcls()
        p = subcmd.add_parser(cmd._alias, help=getattr(cmd, "_help", None))
        p.set_defaults(command=cmd)
        cmd.configure_parser(p)

    return parser
