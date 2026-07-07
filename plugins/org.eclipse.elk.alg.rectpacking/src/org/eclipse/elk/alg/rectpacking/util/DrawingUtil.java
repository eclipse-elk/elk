/*******************************************************************************
 * Copyright (c) 2018, 2020 Kiel University and others.
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
package org.eclipse.elk.alg.rectpacking.util;

import java.util.List;

import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.graph.ElkNode;

/**
 * Class that offers different methods related to coordinates.
 */
public final class DrawingUtil {

    // private constructor
    private DrawingUtil() {
    }

    /**
     * Computes the altered scale measure for a given height and width according to the desired aspect ratio.
     * 
     * @param width The given width.
     * @param height The given height.
     * @param dar The desired aspect ratio.
     * @return altered scale measure calculated the given parameters.
     */
    public static double computeScaleMeasure(final double width, final double height, final double dar) {
        return Math.min(dar / width, 1d / height);
    }

    /**
     * Resets coordinates of every element to (0,0).
     * 
     * @param graph
     *            given graph to reset the coordinates from.
     */
    public static void resetCoordinates(final List<ElkNode> graph) {
        for (ElkNode node : graph) {
            node.setLocation(0, 0);
        }
    }

    /**
     * Calculates the maximum width and height for the given list of {@link RectRow}s.
     * Since the node node spacing is included in the width of a row, it has to be subtracted.
     * @param rows The rows.
     * @param nodeNodeSpacing The spacing between two nodes.
     */
    public static KVector calculateDimensions(final List<RectRow> rows, double nodeNodeSpacing) {
        double maxWidth = 0;
        double newHeight = 0;
        int index = 0;
        for (RectRow row : rows) {
            maxWidth = Math.max(maxWidth, row.getWidth());
            newHeight += row.getHeight() + (index > 0 ? nodeNodeSpacing : 0);
            index++;
        }
        return new KVector(maxWidth, newHeight);
    }

    /**
     * Calculates the maximum width and height for the given list of {@link RectRow}s.
     * Since the node node spacing is included in the width of a row, it has to be subtracted.
     * @param rectangles The rectangles.
     */
    public static KVector calculateDimensions(final List<ElkNode> rectangles) {
        double maxWidth = 0;
        double maxHeight = 0;
        for (ElkNode node : rectangles) {
            maxWidth = Math.max(node.getWidth() + node.getX(), maxWidth);
            maxHeight = Math.max(node.getHeight() + node.getY(), maxHeight);
        }
        return new KVector(maxWidth, maxHeight);
    }
}
