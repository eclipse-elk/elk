/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
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
package org.eclipse.elk.core.ui.views;

import java.util.LinkedList;
import java.util.List;

import org.eclipse.elk.core.service.ILayoutConfigurationStore;
import org.eclipse.elk.core.service.ILayoutSetup;
import org.eclipse.elk.core.service.LayoutConfigurationManager;
import org.eclipse.elk.core.validation.LayoutOptionValidator;
import org.eclipse.ui.IWorkbenchPart;
import org.eclipse.ui.views.properties.IPropertySource;
import org.eclipse.ui.views.properties.IPropertySourceProvider;

import com.google.inject.Inject;
import com.google.inject.Provider;

/**
 * A property source provider used by the layout view.
 * 
 * <p>Subclasses of this class can be bound in an {@link ILayoutSetup} injector for customization.</p>
 *
 * @author msp
 */
public class LayoutPropertySourceProvider implements IPropertySourceProvider {
    
    /**
     * The layout configuration store provider interacting with the created property sources.
     */
    @Inject(optional = true)
    private ILayoutConfigurationStore.Provider configurationStoreProvider;
    
    /**
     * The layout configuration manager used to handle configuration stores.
     */
    @Inject
    private LayoutConfigurationManager configManager;
    
    /**
     * Provider for validators used to check bounds of layout option values.
     */
    @Inject
    private Provider<LayoutOptionValidator> layoutOptionValidatorProvider;

    /**
     * The workbench part containing the current selection.
     */
    private IWorkbenchPart workbenchPart;
    
    /**
     * The requested configuration stores.
     */
    private final List<ILayoutConfigurationStore> configStores = new LinkedList<ILayoutConfigurationStore>();
    
    /**
     * Return the currently tracked workbench part.
     * 
     * @return the current workbench part
     */
    public IWorkbenchPart getWorkbenchPart() {
        return workbenchPart;
    }
    
    /**
     * Set the currently tracked workbench part.
     * 
     * @param theworkbenchPart a workbench part, or {@code null}
     */
    public void setWorkbenchPart(final IWorkbenchPart theworkbenchPart) {
        this.workbenchPart = theworkbenchPart;
    }
    
    /**
     * Return the list of requested layout configuration stores.
     */
    public List<ILayoutConfigurationStore> getConfigurationStores() {
        return configStores;
    }
    
    /**
     * Return the active layout configuration manager.
     */
    public LayoutConfigurationManager getConfigurationManager() {
        return configManager;
    }
    
    @Override
    public IPropertySource getPropertySource(final Object object) {
        if (configurationStoreProvider != null) {
            ILayoutConfigurationStore config = configurationStoreProvider.get(workbenchPart, object);
            if (config != null) {
                configStores.add(config);
                return createPropertySource(config);
            }
        }
        return null;
    }
    
    /**
     * Create a property source for the given configuration store.
     */
    protected IPropertySource createPropertySource(final ILayoutConfigurationStore config) {
        LayoutPropertySource result = new LayoutPropertySource(config, getConfigurationManager());
        result.setValidator(layoutOptionValidatorProvider.get());
        return result;
    }
    
}
