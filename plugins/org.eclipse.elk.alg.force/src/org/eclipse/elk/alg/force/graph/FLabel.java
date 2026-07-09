/*******************************************************************************
 * Copyright (c) 2011, 2020 Kiel University and others.
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
package org.eclipse.elk.alg.force.graph;

import org.eclipse.elk.alg.force.options.ForceOptions;
import org.eclipse.elk.core.math.KVector;

/**
 * A label in the force graph.
 */
public final class FLabel extends FParticle {
    
    /** the serial version UID. */
    private static final long serialVersionUID = 9047772256368142239L;
    
    /** the edge this label is associated to. */
    private FEdge edge;
    /** label text. */
    private String text;
    
    /**
     * Create a new label. The label is also put into the edge's list of labels.
     * 
     * @param text the text of the new label
     * @param fedge edge corresponding to this label
     */
    public FLabel(final FEdge fedge, final String text) {
        this.edge = fedge;
        this.text = text;
        edge.getLabels().add(this);
    }
    
    @Override
    public String toString() {
        if (text == null || text.length() == 0) {
            return "l[" + edge.toString() + "]";
        } else {
            return "l_" + text;
        }
    }
        
    /**
     * Returns the text of this label.
     * 
     * @return text of this label
     */
    public String getText() {
        return text;
    }
    
    /**
     * Returns the associated edge.
     * 
     * @return edge this label is associated to.
     */
    public FEdge getEdge() {
        return edge;
    }
    
    /**
     * Refresh the label position, that is, place it as specified by {@link ForceOptions#EDGE_LABELS_INLINE}.
     */
    public void refreshPosition() {
        
        boolean placeInline = getProperty(ForceOptions.EDGE_LABELS_INLINE);
        KVector src = edge.getSource().getPosition();
        KVector tgt = edge.getTarget().getPosition();

        if (placeInline) {
            // SUPPRESS CHECKSTYLE NEXT 2 MagicNumber
            KVector srcToTgt = tgt.clone().sub(src).scale(0.5);
            KVector toLabelCenter = getSize().clone().scale(0.5);
            KVector newLabelPosition = src.clone().add(srcToTgt).sub(toLabelCenter);
            getPosition().set(newLabelPosition);
        } else {
            // TODO The spacing has to be retrieved from the surrounding node/graph
            double spacing = edge.getProperty(ForceOptions.SPACING_EDGE_LABEL).doubleValue();
            // TODO add support for head and tail labels
            KVector pos = getPosition();
            if (src.x >= tgt.x) {
                if (src.y >= tgt.y) {
                    // CASE1, src top left, tgt bottom right
                    pos.x = tgt.x + ((src.x - tgt.x) / 2) + spacing;
                    pos.y = tgt.y + ((src.y - tgt.y) / 2) - spacing - getSize().y;
                } else {
                    // CASE2, src bottom left, tgt top right
                    pos.x = tgt.x + ((src.x - tgt.x) / 2) + spacing;
                    pos.y = src.y + ((tgt.y - src.y) / 2) + spacing;
                }
            } else {
                if (src.y >= tgt.y) {
                    // CASE2, src top right, tgt bottom left
                    pos.x = src.x + ((tgt.x - src.x) / 2) + spacing;
                    pos.y = tgt.y + ((src.y - tgt.y) / 2) + spacing;
                } else {
                    // CASE1, src bottom right, tgt top left
                    pos.x = src.x + ((tgt.x - src.x) / 2) + spacing;
                    pos.y = src.y + ((tgt.y - src.y) / 2) - spacing - getSize().y;
                }
            }
        }
    }

}
