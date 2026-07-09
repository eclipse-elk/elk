/**
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
 */
package org.eclipse.elk.graph.impl;

import java.util.Map;
import org.eclipse.elk.graph.EMapPropertyHolder;
import org.eclipse.elk.graph.ElkGraphPackage;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.IPropertyHolder;
import org.eclipse.elk.graph.properties.IPropertyValueProxy;
import org.eclipse.emf.common.notify.NotificationChain;
import org.eclipse.emf.common.util.EMap;
import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.EStructuralFeature;
import org.eclipse.emf.ecore.InternalEObject;
import org.eclipse.emf.ecore.impl.MinimalEObjectImpl;
import org.eclipse.emf.ecore.util.EcoreEMap;
import org.eclipse.emf.ecore.util.InternalEList;

/**
 * <!-- begin-user-doc -->
 * An implementation of the model object '<em><b>EMap Property Holder</b></em>'.
 * <!-- end-user-doc -->
 * <p>
 * The following features are implemented:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.graph.impl.EMapPropertyHolderImpl#getProperties <em>Properties</em>}</li>
 * </ul>
 *
 * @generated
 */
public abstract class EMapPropertyHolderImpl extends MinimalEObjectImpl.Container implements EMapPropertyHolder {
    /**
     * The cached value of the '{@link #getProperties() <em>Properties</em>}' map.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #getProperties()
     * @generated
     * @ordered
     */
    protected EMap<IProperty<?>, Object> properties;

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    protected EMapPropertyHolderImpl() {
        super();
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @Override
    protected EClass eStaticClass() {
        return ElkGraphPackage.Literals.EMAP_PROPERTY_HOLDER;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EMap<IProperty<?>, Object> getProperties() {
        if (properties == null) {
            properties = new EcoreEMap<IProperty<?>,Object>(ElkGraphPackage.Literals.ELK_PROPERTY_TO_VALUE_MAP_ENTRY, ElkPropertyToValueMapEntryImpl.class, this, ElkGraphPackage.EMAP_PROPERTY_HOLDER__PROPERTIES);
        }
        return properties;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated NOT
     */
    public <T> IPropertyHolder setProperty(IProperty<? super T> property, T value) {
        if (value == null) {
            getProperties().removeKey(property);
        } else {
            getProperties().put(property, value);
        }
        
        return this;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated NOT
     */
    @SuppressWarnings("unchecked")
    public <T> T getProperty(IProperty<T> property) {
        Object value = getProperties().get(property);
        if (value instanceof IPropertyValueProxy) {
            value = ((IPropertyValueProxy) value).resolveValue(property);
            if (value != null) {
                getProperties().put(property, value);
                return (T) value;
            }
        } else if (value != null) {
            return (T) value;
        }

        // the reason for the side effect below is that if a default value has been returned 
        // and the object is altered by the user, the user expects the altered object to be 
        // the value of the property in case he asks for the property again
        
        T defaultValue = property.getDefault();
        if (defaultValue instanceof Cloneable) {
            // We are now dealing with a clone of the default value which me may safely store away
            // for further modification
            setProperty(property, defaultValue);
        }
        
        return defaultValue;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated NOT
     */
    public boolean hasProperty(IProperty<?> property) {
        return getProperties().containsKey(property);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated NOT
     */
    public IPropertyHolder copyProperties(IPropertyHolder source) {
        if (source == null) {
            return this;
        }

        if (source instanceof EMapPropertyHolder) {
            EMapPropertyHolder other = (EMapPropertyHolder) source;
            EMap<IProperty<?>, Object> ourProps = this.getProperties();
            for (Map.Entry<IProperty<?>, Object> entry : other.getProperties()) {
                Object value = entry.getValue();
                if (value instanceof IPropertyValueProxy) {
                    IPropertyValueProxy proxy = (IPropertyValueProxy) value;
                    Object newValue = proxy.resolveValue(entry.getKey());
                    if (newValue != null) {
                        entry.setValue(newValue);
                        value = newValue;
                    }
                }
                ourProps.put(entry.getKey(), value);
            }
        } else {
            this.getProperties().putAll(source.getAllProperties());
        }

        return this;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated NOT
     */
    public Map<IProperty<?>, Object> getAllProperties() {
        EMap<IProperty<?>, Object> props = getProperties();
        
        // check for unresolved properties
        for (Map.Entry<IProperty<?>, Object> entry : props) {
            if (entry.getValue() instanceof IPropertyValueProxy && entry.getKey() != null) {
                IPropertyValueProxy proxy = (IPropertyValueProxy) entry.getValue();
                
                // Try to resolve the proxy's value, maybe the layout option was 
                // registered by now. If not, we preserve the proxy. 
                Object value = proxy.resolveValue(entry.getKey());
                if (value != null) {
                    entry.setValue(value);
                }
            }
        }
        
        return props.map();
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @Override
    public NotificationChain eInverseRemove(InternalEObject otherEnd, int featureID, NotificationChain msgs) {
        switch (featureID) {
            case ElkGraphPackage.EMAP_PROPERTY_HOLDER__PROPERTIES:
                return ((InternalEList<?>)getProperties()).basicRemove(otherEnd, msgs);
        }
        return super.eInverseRemove(otherEnd, featureID, msgs);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @Override
    public Object eGet(int featureID, boolean resolve, boolean coreType) {
        switch (featureID) {
            case ElkGraphPackage.EMAP_PROPERTY_HOLDER__PROPERTIES:
                if (coreType) return getProperties();
                else return getProperties().map();
        }
        return super.eGet(featureID, resolve, coreType);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @Override
    public void eSet(int featureID, Object newValue) {
        switch (featureID) {
            case ElkGraphPackage.EMAP_PROPERTY_HOLDER__PROPERTIES:
                ((EStructuralFeature.Setting)getProperties()).set(newValue);
                return;
        }
        super.eSet(featureID, newValue);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @Override
    public void eUnset(int featureID) {
        switch (featureID) {
            case ElkGraphPackage.EMAP_PROPERTY_HOLDER__PROPERTIES:
                getProperties().clear();
                return;
        }
        super.eUnset(featureID);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @Override
    public boolean eIsSet(int featureID) {
        switch (featureID) {
            case ElkGraphPackage.EMAP_PROPERTY_HOLDER__PROPERTIES:
                return properties != null && !properties.isEmpty();
        }
        return super.eIsSet(featureID);
    }

} //EMapPropertyHolderImpl
