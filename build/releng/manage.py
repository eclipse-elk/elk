# Copyright (c) 2020 Kiel University and others.
# This program and the accompanying materials are made available under the
# terms of the Eclipse Public License 2.0 which is available at
# http://www.eclipse.org/legal/epl-2.0.
# 
# This Source Code may also be made available under the following Secondary
# Licenses when the conditions for such availability set forth in the Eclipse
# Public License v. 2.0 are satisfied: GPL-3.0 which is available at
# https://www.gnu.org/licenses/gpl-3.0-standalone.html.
#
# SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later

"""Helper script to manage our downloads area.

This script currently only lists files and directories in our downloads area on the
Eclipse download servers. Future versions could support removing files and directories
as well. If someone adds such functionality, it should default to printing what it
would do first, and require an explicit command line switch to make it actually do it.
"""

import argparse
from pathlib import Path
import shutil

# The default directory all path operations will be relative to.
DEFAULT_BASE_DIR = "/home/data/httpd/download.eclipse.org/elk/"




#
#       #  ####  #####
#       # #        #
#       #  ####    #
#       #      #   #
#       # #    #   #
####### #  ####    #

def print_path(path: Path, curr_depth: int, max_depth: int) -> None:
    """Prints the name of the given path to the console, with the correct indentation. If
       the path happens to be a directory and we have not reached the maximum depth yet, we
       call ourself recursively on its content as well."""

    if path.is_file():
        print(("  " * curr_depth) + F"f {path.name}")

    elif path.is_dir():
        # The name of the depth-0 directory should be printed in full
        if curr_depth == 0:
            path_name = str(path)
        else:
            path_name = path.name
        print(("  " * curr_depth) + F"d {path_name}")

        # Check if we need to recurse
        if curr_depth < max_depth:
            for d in [x for x in path.iterdir() if x.is_dir()]:
                print_path(d, curr_depth + 1, max_depth)

            for f in [x for x in path.iterdir() if x.is_file()]:
                print_path(f, curr_depth + 1, max_depth)


def cmd_list(args: argparse.Namespace) -> None:
    print_path(Path(args.base_path), 0, args.depth)




#     #                      
##   ##  ####  #    # ###### 
# # # # #    # #    # #      
#  #  # #    # #    # #####  
#     # #    # #    # #      
#     # #    #  #  #  #      
#     #  ####    ##   ###### 

def cmd_mv(args: argparse.Namespace) -> None:
    shutil.move(args.source, args.dest)




 #####                #                          ######
#     # #    # #####  #       # #    # ######    #     #   ##   #####   ####  # #    #  ####
#       ##  ## #    # #       # ##   # #         #     #  #  #  #    # #      # ##   # #    #
#       # ## # #    # #       # # #  # #####     ######  #    # #    #  ####  # # #  # #
#       #    # #    # #       # #  # # #         #       ###### #####       # # #  # # #  ###
#     # #    # #    # #       # #   ## #         #       #    # #   #  #    # # #   ## #    #
 #####  #    # #####  ####### # #    # ######    #       #    # #    #  ####  # #    #  ####

# Setup the main parser
parser = argparse.ArgumentParser(description="Manages folders on the download server.")
parser.add_argument(
    "--base-path", "-bp",
    default=DEFAULT_BASE_DIR,
    help="Sets the base directory all other commands operate on.")
subparsers = parser.add_subparsers(help="Available commands.")

# The list command
parser_list = subparsers.add_parser(
    "list",
    help="Lists files and folders in and, optionally, below the base folder.")
parser_list.set_defaults(func=cmd_list)
parser_list.add_argument(
    "--depth",
    type=int,
    default=1,
    help="How many folder levels to go down. By default, only the base folder's content is listed.")

# The move command
parser_move = subparsers.add_parser(
    "mv",
    help="Move a folder to another location.")
parser_move.set_defaults(func=cmd_mv)
parser_move.add_argument(
    "source",
    type=str,
    help="Path to the folder to move.")
parser_move.add_argument(
    "dest",
    type=str,
    help="Path to the folder's new parent folder.")

# Go parse stuff!
args = parser.parse_args()
args.func(args)
