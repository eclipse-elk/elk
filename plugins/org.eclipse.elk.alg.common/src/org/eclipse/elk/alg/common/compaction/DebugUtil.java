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
package org.eclipse.elk.alg.common.compaction;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.eclipse.elk.alg.common.compaction.oned.CGraph;
import org.eclipse.elk.alg.common.compaction.oned.CNode;
import org.eclipse.elk.core.math.KVector;

/**
 * Debug utilities for the compaction algorithms.
 */
public final class DebugUtil {

    private DebugUtil() { }
    
    /**
     * For debugging. Writes hitboxes to svg.
     * 
     * @param fileName
     *          filename.
     */
    public static void drawHitboxes(final CGraph cGraph, final String fileName) {

        // determine viewBox
        KVector topLeft = new KVector(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
        KVector bottomRight = new KVector(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
        for (CNode cNode : cGraph.cNodes) {
            topLeft.x = Math.min(topLeft.x, cNode.hitbox.x);
            topLeft.y = Math.min(topLeft.y, cNode.hitbox.y);
            bottomRight.x = Math.max(bottomRight.x, cNode.hitbox.x + cNode.hitbox.width);
            bottomRight.y = Math.max(bottomRight.y, cNode.hitbox.y + cNode.hitbox.height);
        }
        
        KVector size = bottomRight.clone().sub(topLeft);

        // drawing hitboxes to svg
        PrintWriter out;
        try {
            out = new PrintWriter(new FileWriter(fileName));

            out.println("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
            out.println("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"100%\" height=\"100%\""
                    + "  viewBox=\"" + (topLeft.x) + " " + (topLeft.y) 
                    + " " + size.x + " " + size.y + "\">");

            out.println("<defs><marker id=\"markerArrow\" markerWidth=\"10\" "
                    + "markerHeight=\"10\" refX=\"0\" refY=\"3\" orient=\"auto\">"
                    + "  <path d=\"M0,0 L0,6 L9,3 z\" style=\"fill: #000000;\" />"
                    + "</marker></defs>");

            for (CNode cNode : cGraph.cNodes) {
                
                // the node's representation
                out.println(cNode.getDebugSVG());

                // the constraints
                for (CNode inc : cNode.constraints) {
                    out.println("<line x1=\"" + (inc.hitbox.x)
                            + "\" y1=\"" + (inc.hitbox.y + inc.hitbox.height / 2) + "\" x2=\""
                            + (cNode.hitbox.x + cNode.hitbox.width) + "\" y2=\""
                            + (cNode.hitbox.y + cNode.hitbox.height / 2)
                            + "\" stroke=\"black\" opacity=\"0.4\""
                            + " style=\"marker-end: url(#markerArrow);\" />");
                }
            }

            out.println("</svg>");

            out.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
}
