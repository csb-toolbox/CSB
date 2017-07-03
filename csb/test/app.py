"""
CSB Test Runner app. Run with -h to see the app's documentation.
"""

from csb.test import Console


def main():
    return Console('csb.test.cases.*')


if __name__ == '__main__':
    main()
