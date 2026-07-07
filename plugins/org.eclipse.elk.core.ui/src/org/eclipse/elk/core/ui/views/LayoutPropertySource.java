/*******************************************************************************
 * Copyright (c) 2009, 2015 Kiel University and others.
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

import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.eclipse.elk.core.data.LayoutAlgorithmData;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.data.LayoutOptionData;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.service.ILayoutConfigurationStore;
import org.eclipse.elk.core.service.ILayoutSetup;
import org.eclipse.elk.core.service.LayoutConfigurationManager;
import org.eclipse.elk.core.ui.Messages;
import org.eclipse.elk.core.ui.util.ElkUiUtil;
import org.eclipse.elk.core.util.Pair;
import org.eclipse.elk.core.validation.GraphIssue;
import org.eclipse.elk.core.validation.LayoutOptionValidator;
import org.eclipse.ui.views.properties.IPropertyDescriptor;
import org.eclipse.ui.views.properties.IPropertySource;

/**
 * A property source for layout options. This class can be specialized by binding a subclass of
 * {@link LayoutPropertySourceProvider} in an {@link ILayoutSetup} injector.
 *
 * @author msp
 */
public class LayoutPropertySource implements IPropertySource {
    
    /**
     * The layout configuration for this property source.
     */
    private ILayoutConfigurationStore layoutConfig;
    /**
     * Array of property descriptors for the option data.
     */
    private IPropertyDescriptor[] propertyDescriptors;
    /**
     * Set of layout option identifiers that can affect the visibility of other options.
     */
    private final Set<String> dependencyOptions = new HashSet<String>();
    /**
     * The layout configuration manager used to handle configuration stores.
     */
    private final LayoutConfigurationManager configManager;
    /**
     * The validator used to check bounds of layout option values.
     */
    private LayoutOptionValidator layoutOptionValidator;

    /**
     * Creates a layout property source for the given layout configuration.
     */
    public LayoutPropertySource(final ILayoutConfigurationStore config,
            final LayoutConfigurationManager manager) {
        this.layoutConfig = config;
        this.configManager = manager;
    }
    
    /**
     * Set the validator used to check bounds of layout option values.
     */
    public void setValidator(final LayoutOptionValidator validator) {
        this.layoutOptionValidator = validator;
    }
    
    /**
     * Return the layout configuration for this property source.
     */
    protected ILayoutConfigurationStore getLayoutConfig() {
        return layoutConfig;
    }
    
    @Override
    public IPropertyDescriptor[] getPropertyDescriptors() {
        if (propertyDescriptors == null) {
            Set<LayoutOptionData> optionData = configManager.getSupportedOptions(layoutConfig);
            
            // Filter the options hidden by visibility settings and option dependencies
            filterOptions(optionData);
            
            propertyDescriptors = new IPropertyDescriptor[optionData.size()];
            Iterator<LayoutOptionData> optionIter = optionData.iterator();
            for (int i = 0; i < propertyDescriptors.length; ++i) {
                propertyDescriptors[i] = createPropertyDescriptor(optionIter.next());
            }
        }
        return propertyDescriptors;
    }
    
    /**
     * Create a property descriptor for the given option data.
     */
    protected IPropertyDescriptor createPropertyDescriptor(final LayoutOptionData data) {
        return new LayoutPropertyDescriptor(this, data, layoutConfig.getOptionTargets());
    }
    
    /**
     * Remove options that should not be visible from the given list. Options that have dependencies
     * are only visible if the dependencies are met. A dependency is met if it has a target value
     * that equals the actual value, or if it has no target value and the actual value is anything
     * but {@code null}.
     * 
     * <p>Override this method in a subclass in order to further restrict the visible options.</p>
     * 
     * @param optionData a list of option meta data
     */
    protected void filterOptions(final Set<LayoutOptionData> optionData) {
        // the layout algorithm option always affects other options
        dependencyOptions.add(CoreOptions.ALGORITHM.getId());
        
        Iterator<LayoutOptionData> optionIter = optionData.iterator();
        while (optionIter.hasNext()) {
            LayoutOptionData option = optionIter.next();
            boolean visible = option.getVisibility() != LayoutOptionData.Visibility.HIDDEN;
            if (visible) {
                visible = option.getDependencies().isEmpty();
                for (Pair<LayoutOptionData, Object> dependency : option.getDependencies()) {
                    // if at least one dependency is met, the option is made visible
                    LayoutOptionData targetOption = dependency.getFirst();
                    dependencyOptions.add(targetOption.getId());
                    Object expectedValue = dependency.getSecond();
                    Object value = configManager.getOptionValue(targetOption, layoutConfig);
                    if (expectedValue == null && value != null
                            || expectedValue != null && expectedValue.equals(value)) {
                        visible = true;
                        break;
                    }
                }
            }
            
            if (!visible) {
                optionIter.remove();
            }
        }
    }

    @Override
    public Object getPropertyValue(final Object id) {
        LayoutMetaDataService layoutServices = LayoutMetaDataService.getInstance();
        LayoutOptionData optionData = layoutServices.getOptionData((String) id);
        if (optionData != null) {
            Object value = configManager.getOptionValue(optionData, layoutConfig);
            return translateToUI(value, optionData);
        }
        return null;
    }
    
    /**
     * Translate a layout option value into an object that can be handled by the cell
     * editors of the layout view.
     * 
     * @param value a layout option value
     * @param optionData the corresponding layout option data
     * @return a cell editor value
     */
    @SuppressWarnings("rawtypes")
    public Object translateToUI(final Object value, final LayoutOptionData optionData) {
        if (value == null) {
            return "";
        }
        switch (optionData.getType()) {
        case INT:
        case DOUBLE:           // TextCellEditor
            return value.toString();
        case BOOLEAN:         // ComboBoxCellEditor
            if (value instanceof Boolean) {
                return Integer.valueOf(((Boolean) value) ? 1 : 0);
            } else if (value instanceof String) {
                return Integer.valueOf(Boolean.valueOf((String) value) ? 1 : 0);
            } else {
                return value;
            }
        case ENUM:            // ComboBoxCellEditor
            if (value instanceof Enum<?>) {
                return ((Enum<?>) value).ordinal();
            } else if (value instanceof String) {
                String[] choices = optionData.getChoices();
                for (int i = 0; i < choices.length; i++) {
                    if (choices[i].equals(value)) {
                        return i;
                    }
                }
                return 0;
            }
            return value;
        case ENUMSET:         // MultipleOptionsCellEditor
            Set set = (Set) value;
            String[] result = new String[set.size()];
            
            Iterator iterator = set.iterator();
            for (int i = 0; iterator.hasNext(); i++) {
                Object o = iterator.next();
                
                if (o instanceof Enum) {
                    result[i] = ((Enum) o).name();
                } else {
                    result[i] = ((String) o);
                }
            }
            
            return result;
        case OBJECT:          // TextCellEditor
            return value.toString();
        default:
            return value;
        }
    }
    
    @Override
    public void setPropertyValue(final Object id, final Object thevalue) {
        final LayoutOptionData optionData = LayoutMetaDataService.getInstance()
                .getOptionData((String) id);
        if (optionData != null) {
            Runnable modelChange = new Runnable() {
                public void run() {
                    Object value = translateFromUI(thevalue, optionData);
                    layoutConfig.setOptionValue(optionData.getId(), value.toString());
                }
            };
            ElkUiUtil.runModelChange(modelChange, layoutConfig.getEditingDomain(), Messages.getString("elk.ui.11"));
    
            // If the selected option can affect other options, refresh the whole layout view
            if (dependencyOptions.contains(id)) {
                LayoutViewPart layoutView = LayoutViewPart.findView();
                if (layoutView != null) {
                    layoutView.refresh();
                }
            }
        }
    }
    
    /**
     * Translate a layout option value from an object returned by a cell editor.
     * 
     * @param value a cell editor value
     * @param optionData the corresponding layout option data
     * @return a layout option value
     */
    public Object translateFromUI(final Object value, final LayoutOptionData optionData) {
        switch (optionData.getType()) {
        case STRING:
            return (String) value;
        case BOOLEAN:
            return Boolean.valueOf((Integer) value == 1);
        case ENUM:
            return optionData.getEnumValue((Integer) value);
        case ENUMSET:
            // The returned value is a string array that we will turn into a string
            // of elements separated by whitespace. We can then use LayoutOptionData
            // to obtain a proper set
            StringBuilder elementString = new StringBuilder();
            for (String s : (String[]) value) {
                elementString.append(" ").append(s);
            }
            return optionData.parseValue(elementString.toString());
        default:
            return optionData.parseValue((String) value);
        }
    }
    
    /**
     * Validate the given property value and return a list of issues. If the value is ok, the list is empty.
     * If no validator has been configured for this property source, {@code null} is returned.
     */
    public List<GraphIssue> validatePropertyValue(final LayoutOptionData optionData, final Object value) {
        if (layoutOptionValidator == null) {
            return null;
        }
        return layoutOptionValidator.checkProperty(optionData, value, null);
    }
    
    @Override
    public Object getEditableValue() {
        // this feature is currently not required (see interface documentation)
        return null;
    }

    @Override
    public boolean isPropertySet(final Object id) {
        return layoutConfig.getOptionValue((String) id) != null;
    }

    @Override
    public void resetPropertyValue(final Object id) {
        Runnable modelChange = new Runnable() {
            public void run() {
                layoutConfig.setOptionValue((String) id, null);
            }
        };
        ElkUiUtil.runModelChange(modelChange, layoutConfig.getEditingDomain(), Messages.getString("elk.ui.12"));
        
        // if the selected option can affect other options, refresh the whole layout view
        if (dependencyOptions.contains(id)) {
            LayoutViewPart layoutView = LayoutViewPart.findView();
            if (layoutView != null) {
                layoutView.refresh();
            }
        }
    }
    
    /**
     * Returns an identifier for a displayed layout algorithm name. The result is the identifier of
     * an algorithm whose name is a prefix of the displayed name. If there are multiple such
     * algorithms, the one with the longest prefix is taken.
     * 
     * @param displayedName a displayed name of a layout algorithm
     * @return the corresponding identifier, or {@code null} if no match is found
     */
    public static String getLayoutHint(final String displayedName) {
        String bestHint = null;
        int bestLength = 0;
        for (LayoutAlgorithmData layouterData : LayoutMetaDataService.getInstance().getAlgorithmData()) {
            String name = layouterData.getName();
            if (displayedName.startsWith(name) && name.length() > bestLength) {
                bestHint = layouterData.getId();
                bestLength = name.length();
            }
        }
        return bestHint;
    }

}
