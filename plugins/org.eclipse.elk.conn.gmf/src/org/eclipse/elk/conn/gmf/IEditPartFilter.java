/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
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
package org.eclipse.elk.conn.gmf;

import org.eclipse.gef.EditPart;

import com.google.inject.ImplementedBy;

/**
 * Interface for edit part filters. Use this to exclude certain diagram parts from automatic layout.
 */
@ImplementedBy(IEditPartFilter.DefaultImpl.class)
public interface IEditPartFilter {
    
    /**
     * Whether to accept the given edit part and include it in the layout graph.
     */
    boolean filter(EditPart editPart);
    
    /**
     * This implementation includes all edit parts (returns always {@code true}).
     */
    class DefaultImpl implements IEditPartFilter {

        @Override
        public boolean filter(final EditPart editPart) {
            return true;
        }
        
    }

}
