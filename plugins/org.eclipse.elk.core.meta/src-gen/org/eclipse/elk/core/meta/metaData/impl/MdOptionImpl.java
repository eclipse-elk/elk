/**
 * ******************************************************************************
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
 *  ******************************************************************************
 */
package org.eclipse.elk.core.meta.metaData.impl;

import java.util.Collection;

import org.eclipse.elk.core.meta.metaData.MdOption;
import org.eclipse.elk.core.meta.metaData.MdOptionDependency;
import org.eclipse.elk.core.meta.metaData.MdOptionTargetType;
import org.eclipse.elk.core.meta.metaData.MetaDataPackage;

import org.eclipse.emf.common.notify.Notification;
import org.eclipse.emf.common.notify.NotificationChain;

import org.eclipse.emf.common.util.EList;

import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.InternalEObject;

import org.eclipse.emf.ecore.impl.ENotificationImpl;

import org.eclipse.emf.ecore.util.EDataTypeEList;
import org.eclipse.emf.ecore.util.EObjectContainmentEList;
import org.eclipse.emf.ecore.util.InternalEList;

import org.eclipse.xtext.common.types.JvmTypeReference;

import org.eclipse.xtext.xbase.XExpression;

/**
 * <!-- begin-user-doc -->
 * An implementation of the model object '<em><b>Md Option</b></em>'.
 * <!-- end-user-doc -->
 * <p>
 * The following features are implemented:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdOptionImpl#isDeprecated <em>Deprecated</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdOptionImpl#isAdvanced <em>Advanced</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdOptionImpl#isProgrammatic <em>Programmatic</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdOptionImpl#isOutput <em>Output</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdOptionImpl#isGlobal <em>Global</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdOptionImpl#getType <em>Type</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdOptionImpl#getLabel <em>Label</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdOptionImpl#getDescription <em>Description</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdOptionImpl#getDefaultValue <em>Default Value</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdOptionImpl#getLowerBound <em>Lower Bound</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdOptionImpl#getUpperBound <em>Upper Bound</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdOptionImpl#getTargets <em>Targets</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdOptionImpl#getLegacyIds <em>Legacy Ids</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.impl.MdOptionImpl#getDependencies <em>Dependencies</em>}</li>
 * </ul>
 *
 * @generated
 */
public class MdOptionImpl extends MdGroupOrOptionImpl implements MdOption
{
  /**
   * The default value of the '{@link #isDeprecated() <em>Deprecated</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isDeprecated()
   * @generated
   * @ordered
   */
  protected static final boolean DEPRECATED_EDEFAULT = false;

  /**
   * The cached value of the '{@link #isDeprecated() <em>Deprecated</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isDeprecated()
   * @generated
   * @ordered
   */
  protected boolean deprecated = DEPRECATED_EDEFAULT;

  /**
   * The default value of the '{@link #isAdvanced() <em>Advanced</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isAdvanced()
   * @generated
   * @ordered
   */
  protected static final boolean ADVANCED_EDEFAULT = false;

  /**
   * The cached value of the '{@link #isAdvanced() <em>Advanced</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isAdvanced()
   * @generated
   * @ordered
   */
  protected boolean advanced = ADVANCED_EDEFAULT;

  /**
   * The default value of the '{@link #isProgrammatic() <em>Programmatic</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isProgrammatic()
   * @generated
   * @ordered
   */
  protected static final boolean PROGRAMMATIC_EDEFAULT = false;

  /**
   * The cached value of the '{@link #isProgrammatic() <em>Programmatic</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isProgrammatic()
   * @generated
   * @ordered
   */
  protected boolean programmatic = PROGRAMMATIC_EDEFAULT;

  /**
   * The default value of the '{@link #isOutput() <em>Output</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isOutput()
   * @generated
   * @ordered
   */
  protected static final boolean OUTPUT_EDEFAULT = false;

  /**
   * The cached value of the '{@link #isOutput() <em>Output</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isOutput()
   * @generated
   * @ordered
   */
  protected boolean output = OUTPUT_EDEFAULT;

  /**
   * The default value of the '{@link #isGlobal() <em>Global</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isGlobal()
   * @generated
   * @ordered
   */
  protected static final boolean GLOBAL_EDEFAULT = false;

  /**
   * The cached value of the '{@link #isGlobal() <em>Global</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isGlobal()
   * @generated
   * @ordered
   */
  protected boolean global = GLOBAL_EDEFAULT;

  /**
   * The cached value of the '{@link #getType() <em>Type</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getType()
   * @generated
   * @ordered
   */
  protected JvmTypeReference type;

  /**
   * The default value of the '{@link #getLabel() <em>Label</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getLabel()
   * @generated
   * @ordered
   */
  protected static final String LABEL_EDEFAULT = null;

  /**
   * The cached value of the '{@link #getLabel() <em>Label</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getLabel()
   * @generated
   * @ordered
   */
  protected String label = LABEL_EDEFAULT;

  /**
   * The default value of the '{@link #getDescription() <em>Description</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getDescription()
   * @generated
   * @ordered
   */
  protected static final String DESCRIPTION_EDEFAULT = null;

  /**
   * The cached value of the '{@link #getDescription() <em>Description</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getDescription()
   * @generated
   * @ordered
   */
  protected String description = DESCRIPTION_EDEFAULT;

  /**
   * The cached value of the '{@link #getDefaultValue() <em>Default Value</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getDefaultValue()
   * @generated
   * @ordered
   */
  protected XExpression defaultValue;

  /**
   * The cached value of the '{@link #getLowerBound() <em>Lower Bound</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getLowerBound()
   * @generated
   * @ordered
   */
  protected XExpression lowerBound;

  /**
   * The cached value of the '{@link #getUpperBound() <em>Upper Bound</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getUpperBound()
   * @generated
   * @ordered
   */
  protected XExpression upperBound;

  /**
   * The cached value of the '{@link #getTargets() <em>Targets</em>}' attribute list.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getTargets()
   * @generated
   * @ordered
   */
  protected EList<MdOptionTargetType> targets;

  /**
   * The cached value of the '{@link #getLegacyIds() <em>Legacy Ids</em>}' attribute list.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getLegacyIds()
   * @generated
   * @ordered
   */
  protected EList<String> legacyIds;

  /**
   * The cached value of the '{@link #getDependencies() <em>Dependencies</em>}' containment reference list.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getDependencies()
   * @generated
   * @ordered
   */
  protected EList<MdOptionDependency> dependencies;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  protected MdOptionImpl()
  {
    super();
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  protected EClass eStaticClass()
  {
    return MetaDataPackage.Literals.MD_OPTION;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public boolean isDeprecated()
  {
    return deprecated;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setDeprecated(boolean newDeprecated)
  {
    boolean oldDeprecated = deprecated;
    deprecated = newDeprecated;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_OPTION__DEPRECATED, oldDeprecated, deprecated));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public boolean isAdvanced()
  {
    return advanced;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setAdvanced(boolean newAdvanced)
  {
    boolean oldAdvanced = advanced;
    advanced = newAdvanced;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_OPTION__ADVANCED, oldAdvanced, advanced));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public boolean isProgrammatic()
  {
    return programmatic;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setProgrammatic(boolean newProgrammatic)
  {
    boolean oldProgrammatic = programmatic;
    programmatic = newProgrammatic;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_OPTION__PROGRAMMATIC, oldProgrammatic, programmatic));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public boolean isOutput()
  {
    return output;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setOutput(boolean newOutput)
  {
    boolean oldOutput = output;
    output = newOutput;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_OPTION__OUTPUT, oldOutput, output));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public boolean isGlobal()
  {
    return global;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setGlobal(boolean newGlobal)
  {
    boolean oldGlobal = global;
    global = newGlobal;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_OPTION__GLOBAL, oldGlobal, global));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public JvmTypeReference getType()
  {
    return type;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public NotificationChain basicSetType(JvmTypeReference newType, NotificationChain msgs)
  {
    JvmTypeReference oldType = type;
    type = newType;
    if (eNotificationRequired())
    {
      ENotificationImpl notification = new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_OPTION__TYPE, oldType, newType);
      if (msgs == null) msgs = notification; else msgs.add(notification);
    }
    return msgs;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setType(JvmTypeReference newType)
  {
    if (newType != type)
    {
      NotificationChain msgs = null;
      if (type != null)
        msgs = ((InternalEObject)type).eInverseRemove(this, EOPPOSITE_FEATURE_BASE - MetaDataPackage.MD_OPTION__TYPE, null, msgs);
      if (newType != null)
        msgs = ((InternalEObject)newType).eInverseAdd(this, EOPPOSITE_FEATURE_BASE - MetaDataPackage.MD_OPTION__TYPE, null, msgs);
      msgs = basicSetType(newType, msgs);
      if (msgs != null) msgs.dispatch();
    }
    else if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_OPTION__TYPE, newType, newType));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public String getLabel()
  {
    return label;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setLabel(String newLabel)
  {
    String oldLabel = label;
    label = newLabel;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_OPTION__LABEL, oldLabel, label));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public String getDescription()
  {
    return description;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setDescription(String newDescription)
  {
    String oldDescription = description;
    description = newDescription;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_OPTION__DESCRIPTION, oldDescription, description));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public XExpression getDefaultValue()
  {
    return defaultValue;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public NotificationChain basicSetDefaultValue(XExpression newDefaultValue, NotificationChain msgs)
  {
    XExpression oldDefaultValue = defaultValue;
    defaultValue = newDefaultValue;
    if (eNotificationRequired())
    {
      ENotificationImpl notification = new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_OPTION__DEFAULT_VALUE, oldDefaultValue, newDefaultValue);
      if (msgs == null) msgs = notification; else msgs.add(notification);
    }
    return msgs;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setDefaultValue(XExpression newDefaultValue)
  {
    if (newDefaultValue != defaultValue)
    {
      NotificationChain msgs = null;
      if (defaultValue != null)
        msgs = ((InternalEObject)defaultValue).eInverseRemove(this, EOPPOSITE_FEATURE_BASE - MetaDataPackage.MD_OPTION__DEFAULT_VALUE, null, msgs);
      if (newDefaultValue != null)
        msgs = ((InternalEObject)newDefaultValue).eInverseAdd(this, EOPPOSITE_FEATURE_BASE - MetaDataPackage.MD_OPTION__DEFAULT_VALUE, null, msgs);
      msgs = basicSetDefaultValue(newDefaultValue, msgs);
      if (msgs != null) msgs.dispatch();
    }
    else if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_OPTION__DEFAULT_VALUE, newDefaultValue, newDefaultValue));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public XExpression getLowerBound()
  {
    return lowerBound;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public NotificationChain basicSetLowerBound(XExpression newLowerBound, NotificationChain msgs)
  {
    XExpression oldLowerBound = lowerBound;
    lowerBound = newLowerBound;
    if (eNotificationRequired())
    {
      ENotificationImpl notification = new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_OPTION__LOWER_BOUND, oldLowerBound, newLowerBound);
      if (msgs == null) msgs = notification; else msgs.add(notification);
    }
    return msgs;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setLowerBound(XExpression newLowerBound)
  {
    if (newLowerBound != lowerBound)
    {
      NotificationChain msgs = null;
      if (lowerBound != null)
        msgs = ((InternalEObject)lowerBound).eInverseRemove(this, EOPPOSITE_FEATURE_BASE - MetaDataPackage.MD_OPTION__LOWER_BOUND, null, msgs);
      if (newLowerBound != null)
        msgs = ((InternalEObject)newLowerBound).eInverseAdd(this, EOPPOSITE_FEATURE_BASE - MetaDataPackage.MD_OPTION__LOWER_BOUND, null, msgs);
      msgs = basicSetLowerBound(newLowerBound, msgs);
      if (msgs != null) msgs.dispatch();
    }
    else if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_OPTION__LOWER_BOUND, newLowerBound, newLowerBound));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public XExpression getUpperBound()
  {
    return upperBound;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public NotificationChain basicSetUpperBound(XExpression newUpperBound, NotificationChain msgs)
  {
    XExpression oldUpperBound = upperBound;
    upperBound = newUpperBound;
    if (eNotificationRequired())
    {
      ENotificationImpl notification = new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_OPTION__UPPER_BOUND, oldUpperBound, newUpperBound);
      if (msgs == null) msgs = notification; else msgs.add(notification);
    }
    return msgs;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setUpperBound(XExpression newUpperBound)
  {
    if (newUpperBound != upperBound)
    {
      NotificationChain msgs = null;
      if (upperBound != null)
        msgs = ((InternalEObject)upperBound).eInverseRemove(this, EOPPOSITE_FEATURE_BASE - MetaDataPackage.MD_OPTION__UPPER_BOUND, null, msgs);
      if (newUpperBound != null)
        msgs = ((InternalEObject)newUpperBound).eInverseAdd(this, EOPPOSITE_FEATURE_BASE - MetaDataPackage.MD_OPTION__UPPER_BOUND, null, msgs);
      msgs = basicSetUpperBound(newUpperBound, msgs);
      if (msgs != null) msgs.dispatch();
    }
    else if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, MetaDataPackage.MD_OPTION__UPPER_BOUND, newUpperBound, newUpperBound));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EList<MdOptionTargetType> getTargets()
  {
    if (targets == null)
    {
      targets = new EDataTypeEList<MdOptionTargetType>(MdOptionTargetType.class, this, MetaDataPackage.MD_OPTION__TARGETS);
    }
    return targets;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EList<String> getLegacyIds()
  {
    if (legacyIds == null)
    {
      legacyIds = new EDataTypeEList<String>(String.class, this, MetaDataPackage.MD_OPTION__LEGACY_IDS);
    }
    return legacyIds;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public EList<MdOptionDependency> getDependencies()
  {
    if (dependencies == null)
    {
      dependencies = new EObjectContainmentEList<MdOptionDependency>(MdOptionDependency.class, this, MetaDataPackage.MD_OPTION__DEPENDENCIES);
    }
    return dependencies;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public NotificationChain eInverseRemove(InternalEObject otherEnd, int featureID, NotificationChain msgs)
  {
    switch (featureID)
    {
      case MetaDataPackage.MD_OPTION__TYPE:
        return basicSetType(null, msgs);
      case MetaDataPackage.MD_OPTION__DEFAULT_VALUE:
        return basicSetDefaultValue(null, msgs);
      case MetaDataPackage.MD_OPTION__LOWER_BOUND:
        return basicSetLowerBound(null, msgs);
      case MetaDataPackage.MD_OPTION__UPPER_BOUND:
        return basicSetUpperBound(null, msgs);
      case MetaDataPackage.MD_OPTION__DEPENDENCIES:
        return ((InternalEList<?>)getDependencies()).basicRemove(otherEnd, msgs);
    }
    return super.eInverseRemove(otherEnd, featureID, msgs);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public Object eGet(int featureID, boolean resolve, boolean coreType)
  {
    switch (featureID)
    {
      case MetaDataPackage.MD_OPTION__DEPRECATED:
        return isDeprecated();
      case MetaDataPackage.MD_OPTION__ADVANCED:
        return isAdvanced();
      case MetaDataPackage.MD_OPTION__PROGRAMMATIC:
        return isProgrammatic();
      case MetaDataPackage.MD_OPTION__OUTPUT:
        return isOutput();
      case MetaDataPackage.MD_OPTION__GLOBAL:
        return isGlobal();
      case MetaDataPackage.MD_OPTION__TYPE:
        return getType();
      case MetaDataPackage.MD_OPTION__LABEL:
        return getLabel();
      case MetaDataPackage.MD_OPTION__DESCRIPTION:
        return getDescription();
      case MetaDataPackage.MD_OPTION__DEFAULT_VALUE:
        return getDefaultValue();
      case MetaDataPackage.MD_OPTION__LOWER_BOUND:
        return getLowerBound();
      case MetaDataPackage.MD_OPTION__UPPER_BOUND:
        return getUpperBound();
      case MetaDataPackage.MD_OPTION__TARGETS:
        return getTargets();
      case MetaDataPackage.MD_OPTION__LEGACY_IDS:
        return getLegacyIds();
      case MetaDataPackage.MD_OPTION__DEPENDENCIES:
        return getDependencies();
    }
    return super.eGet(featureID, resolve, coreType);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @SuppressWarnings("unchecked")
  @Override
  public void eSet(int featureID, Object newValue)
  {
    switch (featureID)
    {
      case MetaDataPackage.MD_OPTION__DEPRECATED:
        setDeprecated((Boolean)newValue);
        return;
      case MetaDataPackage.MD_OPTION__ADVANCED:
        setAdvanced((Boolean)newValue);
        return;
      case MetaDataPackage.MD_OPTION__PROGRAMMATIC:
        setProgrammatic((Boolean)newValue);
        return;
      case MetaDataPackage.MD_OPTION__OUTPUT:
        setOutput((Boolean)newValue);
        return;
      case MetaDataPackage.MD_OPTION__GLOBAL:
        setGlobal((Boolean)newValue);
        return;
      case MetaDataPackage.MD_OPTION__TYPE:
        setType((JvmTypeReference)newValue);
        return;
      case MetaDataPackage.MD_OPTION__LABEL:
        setLabel((String)newValue);
        return;
      case MetaDataPackage.MD_OPTION__DESCRIPTION:
        setDescription((String)newValue);
        return;
      case MetaDataPackage.MD_OPTION__DEFAULT_VALUE:
        setDefaultValue((XExpression)newValue);
        return;
      case MetaDataPackage.MD_OPTION__LOWER_BOUND:
        setLowerBound((XExpression)newValue);
        return;
      case MetaDataPackage.MD_OPTION__UPPER_BOUND:
        setUpperBound((XExpression)newValue);
        return;
      case MetaDataPackage.MD_OPTION__TARGETS:
        getTargets().clear();
        getTargets().addAll((Collection<? extends MdOptionTargetType>)newValue);
        return;
      case MetaDataPackage.MD_OPTION__LEGACY_IDS:
        getLegacyIds().clear();
        getLegacyIds().addAll((Collection<? extends String>)newValue);
        return;
      case MetaDataPackage.MD_OPTION__DEPENDENCIES:
        getDependencies().clear();
        getDependencies().addAll((Collection<? extends MdOptionDependency>)newValue);
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
  public void eUnset(int featureID)
  {
    switch (featureID)
    {
      case MetaDataPackage.MD_OPTION__DEPRECATED:
        setDeprecated(DEPRECATED_EDEFAULT);
        return;
      case MetaDataPackage.MD_OPTION__ADVANCED:
        setAdvanced(ADVANCED_EDEFAULT);
        return;
      case MetaDataPackage.MD_OPTION__PROGRAMMATIC:
        setProgrammatic(PROGRAMMATIC_EDEFAULT);
        return;
      case MetaDataPackage.MD_OPTION__OUTPUT:
        setOutput(OUTPUT_EDEFAULT);
        return;
      case MetaDataPackage.MD_OPTION__GLOBAL:
        setGlobal(GLOBAL_EDEFAULT);
        return;
      case MetaDataPackage.MD_OPTION__TYPE:
        setType((JvmTypeReference)null);
        return;
      case MetaDataPackage.MD_OPTION__LABEL:
        setLabel(LABEL_EDEFAULT);
        return;
      case MetaDataPackage.MD_OPTION__DESCRIPTION:
        setDescription(DESCRIPTION_EDEFAULT);
        return;
      case MetaDataPackage.MD_OPTION__DEFAULT_VALUE:
        setDefaultValue((XExpression)null);
        return;
      case MetaDataPackage.MD_OPTION__LOWER_BOUND:
        setLowerBound((XExpression)null);
        return;
      case MetaDataPackage.MD_OPTION__UPPER_BOUND:
        setUpperBound((XExpression)null);
        return;
      case MetaDataPackage.MD_OPTION__TARGETS:
        getTargets().clear();
        return;
      case MetaDataPackage.MD_OPTION__LEGACY_IDS:
        getLegacyIds().clear();
        return;
      case MetaDataPackage.MD_OPTION__DEPENDENCIES:
        getDependencies().clear();
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
  public boolean eIsSet(int featureID)
  {
    switch (featureID)
    {
      case MetaDataPackage.MD_OPTION__DEPRECATED:
        return deprecated != DEPRECATED_EDEFAULT;
      case MetaDataPackage.MD_OPTION__ADVANCED:
        return advanced != ADVANCED_EDEFAULT;
      case MetaDataPackage.MD_OPTION__PROGRAMMATIC:
        return programmatic != PROGRAMMATIC_EDEFAULT;
      case MetaDataPackage.MD_OPTION__OUTPUT:
        return output != OUTPUT_EDEFAULT;
      case MetaDataPackage.MD_OPTION__GLOBAL:
        return global != GLOBAL_EDEFAULT;
      case MetaDataPackage.MD_OPTION__TYPE:
        return type != null;
      case MetaDataPackage.MD_OPTION__LABEL:
        return LABEL_EDEFAULT == null ? label != null : !LABEL_EDEFAULT.equals(label);
      case MetaDataPackage.MD_OPTION__DESCRIPTION:
        return DESCRIPTION_EDEFAULT == null ? description != null : !DESCRIPTION_EDEFAULT.equals(description);
      case MetaDataPackage.MD_OPTION__DEFAULT_VALUE:
        return defaultValue != null;
      case MetaDataPackage.MD_OPTION__LOWER_BOUND:
        return lowerBound != null;
      case MetaDataPackage.MD_OPTION__UPPER_BOUND:
        return upperBound != null;
      case MetaDataPackage.MD_OPTION__TARGETS:
        return targets != null && !targets.isEmpty();
      case MetaDataPackage.MD_OPTION__LEGACY_IDS:
        return legacyIds != null && !legacyIds.isEmpty();
      case MetaDataPackage.MD_OPTION__DEPENDENCIES:
        return dependencies != null && !dependencies.isEmpty();
    }
    return super.eIsSet(featureID);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public String toString()
  {
    if (eIsProxy()) return super.toString();

    StringBuilder result = new StringBuilder(super.toString());
    result.append(" (deprecated: ");
    result.append(deprecated);
    result.append(", advanced: ");
    result.append(advanced);
    result.append(", programmatic: ");
    result.append(programmatic);
    result.append(", output: ");
    result.append(output);
    result.append(", global: ");
    result.append(global);
    result.append(", label: ");
    result.append(label);
    result.append(", description: ");
    result.append(description);
    result.append(", targets: ");
    result.append(targets);
    result.append(", legacyIds: ");
    result.append(legacyIds);
    result.append(')');
    return result.toString();
  }

} //MdOptionImpl
