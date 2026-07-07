/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
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
package org.eclipse.elk.core.service;

import java.util.Collection;
import java.util.Set;

import org.eclipse.elk.core.data.LayoutOptionData;
import org.eclipse.emf.edit.domain.EditingDomain;
import org.eclipse.ui.IWorkbenchPart;

/**
 * Persistence layer for layout configuration. Implementations are used by the Layout view to modify
 * values of layout options for selected diagram elements. An instance of such an implementation should
 * always be linked to a specific context, e.g. a selected node of a diagram. The {@link #getParent()}
 * method allows to navigate to the parent context, so it should return a new instance (or {@code null}
 * if this is already the top-level context).
 * 
 * <p>Implementations of this interface are not bound directly in an {@link ILayoutSetup} injector
 * because each instance is context-sensitive. Use an implementation of {@link Provider} to include
 * your configuration store in the layout setup.</p>
 *
 * @author msp
 */
public interface ILayoutConfigurationStore {
    
    /**
     * Return the editing domain for applying changes to the stored option values, or {@code null}
     * if no editing domain is necessary.
     */
    EditingDomain getEditingDomain();
    
    /**
     * Get the current value for a layout option in the associated context. The return value can be
     * either a literal value or a serialized form as a String.
     * 
     * @param optionId a layout option identifier
     * @return the layout option value, or {@code null} if no value is associated with the context
     */
    Object getOptionValue(String optionId);
    
    /**
     * Set a new value for a layout option in the associated context.
     * 
     * @param optionId a layout option identifier
     * @param value the new layout option value, or {@code null} if the current value shall be removed
     */
    void setOptionValue(String optionId, String value);
    
    /**
     * Retrieve all options that are affected by this configuration store.
     */
    Collection<String> getAffectedOptions();
    
    /**
     * Determine the applicable layout option targets for the associated context.
     * 
     * @return the applicable layout option targets
     */
    Set<LayoutOptionData.Target> getOptionTargets();
    
    /**
     * Return a configuration store for the parent element of the associated context, if any.
     * This should result in a container element that includes the current context.
     * 
     * @return a configuration store for the parent element or {@code null}
     */
    ILayoutConfigurationStore getParent();
    
    /**
     * Interface for providers of layout configuration stores. Implementations can be bound in an
     * {@link ILayoutSetup} injector for linking them to a {@link IDiagramLayoutConnector} implementation.
     */
    interface Provider {
        
        /**
         * Return a layout configuration store that is able to read and write layout options
         * through annotations of the diagram. This configurator is necessary for the Layout View.
         * 
         * @param workbenchPart a workbench part, or {@code null}
         * @param context a context for layout configuration, usually a selected diagram element
         * @return a layout configuration store, or {@code null} if the given context is not applicable
         */
        ILayoutConfigurationStore get(IWorkbenchPart workbenchPart, Object context);
        
    }

}
