/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.common.polyomino;

import org.eclipse.elk.alg.common.polyomino.structures.Polyomino;

/**
 * Utility class providing a simple filling algorithm for polyominoes with holes whose cells are at least 8-connected
 * (otherwise application of the filling method might result in the filling of unintended cells).<br>
 * Examples<br>
 * <br>
 * 8-connected polyomino:<br>
 * {@code _X_}<br>
 * {@code XOX}<br>
 * {@code _X_},<br>
 * where {@code _} represents an empty cell, {@code X} a filled polyomino cell, and {@code O} a newly filled cell after
 * applying the algorithm.<br>
 * <br>
 * unconnected polyomino:<br>
 * {@code _XXX_X}<br>
 * {@code XOX!_X}<br>
 * {@code _XXX__},<br>
 * where the additional {@code !} represents an accidentally newly filled cell which isn't part of a hole, actually.
 */
public final class ProfileFill {

    /**
     * Constructor should never be called as this is simply a utility class.
     */
    private ProfileFill() {
    }

    /**
     * Fills the holes of a polyomino.
     * @param poly polyomino to be filled
     */
    public static void fillPolyomino(final Polyomino poly) {
        int width = poly.getWidth();
        int[] northProfile = new int[width];
        int[] southProfile = new int[width];

        int height = poly.getHeight();
        int[] eastProfile = new int[height];
        int[] westProfile = new int[height];

        // Compute visible profiles from all four cardinal directions.
        for (int xi = 0; xi < width; xi++) {
            int y = 0;
            while (y < height && !poly.isBlocked(xi, y)) {
                y++;
            }
            northProfile[xi] = y;
        }

        for (int xi = 0; xi < width; xi++) {
            int y = height - 1;
            while (y >= 0 && !poly.isBlocked(xi, y)) {
                y--;
            }
            southProfile[xi] = y;
        }

        for (int yi = 0; yi < height; yi++) {
            int x = 0;
            while (x < width && !poly.isBlocked(x, yi)) {
                x++;
            }
            eastProfile[yi] = x;
        }

        for (int yi = 0; yi < height; yi++) {
            int x = width - 1;
            while (x >= 0 && !poly.isBlocked(x, yi)) {
                x--;
            }
            westProfile[yi] = x;
        }

        // Fill all cells between these computed bounds.
        for (int xi = 0; xi < width; xi++) {
            for (int yi = 0; yi < height; yi++) {
                if (xi < westProfile[yi] && xi > eastProfile[yi] && yi < southProfile[xi] && yi > northProfile[xi]) {
                    poly.setBlocked(xi, yi);
                }
            }
        }
    }

}
