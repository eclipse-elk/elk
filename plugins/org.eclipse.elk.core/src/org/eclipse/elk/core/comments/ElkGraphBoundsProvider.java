/*******************************************************************************
 * Copyright (c) 2016, 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.comments;

import java.awt.geom.Rectangle2D;
import java.awt.geom.Rectangle2D.Double;

import org.eclipse.elk.graph.ElkShape;

/**
 * An {@link IBoundsProvider} that simply returns the bounds defined in an ELK node.
 */
public class ElkGraphBoundsProvider implements IBoundsProvider<ElkShape, ElkShape> {

    @Override
    public Double boundsForTarget(final ElkShape node) {
        return new Rectangle2D.Double(node.getX(), node.getY(), node.getWidth(), node.getHeight());
    }
    
    @Override
    public Double boundsForComment(final ElkShape node) {
        return new Rectangle2D.Double(node.getX(), node.getY(), node.getWidth(), node.getHeight());
    }

}
