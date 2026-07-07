/*******************************************************************************
 * Copyright (c) 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.labels;

import org.eclipse.elk.core.math.KVector;

/**
 * Knows how to resize a label to a given target width.
 * 
 * <p>
 * During layout, a layout algorithm may discover that the size of a label considerably increases
 * the size of a diagram. If the label (or a parent) has a label manager attached, the manager may
 * be called to try and shorten the label to a given target width. The new size is returned to the
 * layout algorithm to work with. Of course, the changes to the label's text need to actually be
 * applied after automatic layout. How that works depends on the visualization framework used.
 * </p>
 */
public interface ILabelManager {
    
    /**
     * Tries to shorten the label to keep it narrower than the given target width. This may increase
     * the label's height. The new dimensions of the label are returned for the layout algorithm to
     * work with.
     * 
     * @param label
     *            the label to shorten.
     * @param targetWidth
     *            the width the label's new dimensions should try not to exceed.
     * @return the label's dimensions after shortening or {@code null} if the label has not been
     *         shortened.
     */
    KVector manageLabelSize(Object label, double targetWidth);
    
}
