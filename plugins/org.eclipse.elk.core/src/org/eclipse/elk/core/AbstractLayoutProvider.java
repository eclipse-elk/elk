/*******************************************************************************
 * Copyright (c) 2008, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core;

/**
 * A layout provider executes a layout algorithm to layout the child elements of a node.
 * <p>
 * 
 * @author ars
 * @author msp
 */
public abstract class AbstractLayoutProvider implements IGraphLayoutEngine {

    /**
     * Initialize the layout provider with the given parameter.
     * 
     * @param parameter
     *            a string used to parameterize the layout provider instance. Most layout providers will have no use
     *            for this parameter. However, some may use it to change their behavior. For example, the GraphViz
     *            library provides a lot of layout algorithms, but we only have a single layout provider of which we
     *            create one instance for each layout algorithm GraphViz provides. This parameter is then used to
     *            control which layout algorithm the layout provider provides access to.
     */
    public void initialize(final String parameter) {
        // do nothing - override in subclasses
    }

    /**
     * Dispose the layout provider by releasing any resources that are held.
     */
    public void dispose() {
        // do nothing - override in subclasses
    }

}
